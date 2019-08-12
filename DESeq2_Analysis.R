#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# RNASEQ Analysis to accompany:
# Thyroid deficiency in utero stimulates white adipose tissue growth and impairs brown adipose thermogenic capacity
#
# Forhead et al, 2019
#
# Sheep
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/2019_Alison_J_Forhead
#
# CTR Code: CTR_ajf1005_0001
#
# RNASeq Analysis Performed by Xiaohui Zhao
# CTR Bioinformatics
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#------------------------------------------------------------------------------

message("+--- Loading in the libraries (start up messages supressed) ---+")
suppressPackageStartupMessages({
  rm(list=ls())
  library('DESeq2')
  library('ggplot2')
  library('RColorBrewer')
  library('gplots')
  library("cowplot")
  library("pheatmap")
  library("ggrepel")
  library("reshape2")
  library("biomaRt")
  library("matrixStats")
  library("plyr")
  library("dplyr")
  library("calibrate")
  library('ggplot2')
  library('ggthemes')
  library('ggalt')
  library('UpSetR')
  library('ComplexHeatmap')
 })


message("+------  Set up some constants e.g. base directories   -------+")

Project  <- "CTR_ajf1005_0001"
Base.dir <- "/Users/xz289/Documents/CTR_ajf1005_0001"
FTSeq.dir <- paste0(Base.dir,"/FTSEQ-COUNTS")
out.dir <- paste0(Base.dir, "/Figures_Tables")
script.dir <- paste0(Base.dir,"/Scripts")
med.dir <- paste0(Base.dir, "/Figures_Tables_2")

elementTextSize <- 10
significance <- 0.05
l2fc <- 1
#source(paste0(script.dir,"/Utility_function.R"))
source(paste0(script.dir, "/Backstage_Functions_19_07_2019.R"))
message("+------------     Set up the sample table    -----------------------------------+")

sampleTable <- read.csv(paste0(Base.dir, "/SampleTable.csv"), header = T)
sampleFiles      <- grep('*_fcounts.txt',list.files(FTSeq.dir),value=TRUE)

## match the sample files with the Data_Name
sampleNames.0      <- gsub("_1_trimmed_Oar_v3.1.star.bam.txt", "", sampleFiles)
sampleNames.1      <- unlist(lapply(sampleNames.0, function(x) strsplit(x, split= "[.]")[[1]][2]))
sampleNames.len <- length(sampleNames.1)
select.seq <- seq(1, sampleNames.len, by = 2)
sampleNames <- sampleNames.1[select.seq]
sampleFilesNames <- cbind(sampleFiles, sampleNames)
colnames(sampleFilesNames) <- c("sampleFiles","Data_Name")

sampleTable.1  <- merge(sampleTable, sampleFilesNames, by="Data_Name")

## refine the sample Table
sampleBarcodes <- paste0("s_", sampleTable.1$Sample.ID,"_", rep(c(3,4), length=74))
sampleType <- sampleTable.1$Treatment
samplePhenotype  <- as.character(sampleTable.1$Gestational.age)
sampleSex <- as.character(sampleTable.1$sex_of_fetus)
Labels <-  paste0("s_", sampleTable.1$Sample.ID)

sampleTable.new <- data.frame(sampleName=sampleBarcodes, fileName=sampleFiles,
                              condition=sampleType, GestationalAge=samplePhenotype,
                              sex=sampleSex, Labels=Labels)
print(sampleTable.new)
write.table(sampleTable.new, file =  paste0(out.dir, "/DEseq2_sample_Table_noDrop.txt"),
            row.names=F, col.names=T)


message("+------------  Removing samples identified by QC  -------------+")

sampleTable.rm <- subset(sampleTable.new, Labels!= "s_3" & Labels != "s_27")
print(sampleTable.rm)
sampleTable.rm$CombineCondition <- paste0(sampleTable.rm$condition, "_", sampleTable.rm$GestationalAge)
message("+------------       Retrieve ensEMBL annotations   ------------+")

ensembl    =  useEnsembl(biomart="ensembl", dataset="oaries_gene_ensembl")
## oaries_gene_ensembl; Sheep genes (Oar_v3.1); Oar_v3.1
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl)
head(ensEMBL2id)

message("------------     Creat ddsHTSeq/dds object       ---------------+")

DESeq2Version <- "DESeq2_1.22.2"

ddsHTSeq.rm <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable.rm, directory=FTSeq.dir, design=~0+condition, betaPrior = F)
ddsHTSeq.rm

dds.rm <- DESeq(ddsHTSeq.rm)
dds.rm <- estimateSizeFactors(dds.rm)
rld.rm <- rlogTransformation(dds.rm, blind=F)
vsd.rm <- varianceStabilizingTransformation(dds.rm, blind=F)

message("+--------------------------      Create results object                                      ------------------+")

normCounts.rm.ori <- counts(dds.rm, normalized=TRUE)
normCounts.rm <- cbind(normCounts.rm.ori, rownames(normCounts.rm.ori))
colnames(normCounts.rm) <- c(colnames(normCounts.rm.ori), "ensembl_gene_id")

normCounts.rm.annot <- merge(normCounts.rm,ensEMBL2id, by="ensembl_gene_id" )

save(dds.rm, rld.rm, vsd.rm, normCounts.rm.annot, file = paste0(out.dir, "/", Project, "-DESeq2_rm3_rm27_dds_rld_normCannot.RData"))

message("+---- Collapse the replicates and creat sample table, dds and results object  ----------+")

sampleTable.collapsed.rm            <- sampleTable.rm %>% group_by(Labels) %>% summarise_all(funs(paste(unique(.), collapse=",")))
sampleTable.collapsed.rm$sampleName <- paste0(sampleTable.collapsed.rm$Labels, "_",sampleTable.collapsed.rm$condition,"_",
                                              sampleTable.collapsed.rm$GestationalAge)
sampleTable.collapsed.rm$CombineCondition <- paste0(sampleTable.collapsed.rm$condition, "_",sampleTable.collapsed.rm$GestationalAge)
write.csv(sampleTable.collapsed.rm, file = paste0(out.dir, "/", Project, "-SampleTable_Collapsed_rm27_rm3.csv"), row.names=F)

dds.collapsed.rm <- collapseReplicates(ddsHTSeq.rm, groupby= ddsHTSeq.rm$Labels,  renameCols = TRUE)
dds.collapsed.rm <- DESeq(dds.collapsed.rm)
dds.collapsed.rm <- estimateSizeFactors(dds.collapsed.rm)

resultsNames(dds.collapsed.rm)
colData(dds.collapsed.rm)
colnames(dds.collapsed.rm)

res_TX_SHAM <- results(dds.collapsed.rm, contrast=c("condition", "TX", "SHAM"))
res_TX_SHAM.df                 <- as.data.frame(res_TX_SHAM)
res_TX_SHAM.df$ensembl_gene_id <- rownames(res_TX_SHAM.df)

res_TX_SHAM.ann                <- merge(res_TX_SHAM.df,ensEMBL2id, by="ensembl_gene_id")
res_TX_SHAM.ann$description    <- gsub("..Source.*", "", res_TX_SHAM.ann$description)
res_TX_SHAM.ann <- subset(res_TX_SHAM.ann, !is.na(log2FoldChange))
## 25056 * 10

res_TX_SHAM.ann.sig            <- subset(res_TX_SHAM.ann, padj <= significance  & abs(log2FoldChange) >= l2fc)
head(res_TX_SHAM.ann.sig)
sigNum <- dim(res_TX_SHAM.ann.sig)[1]
## 1472

save(res_TX_SHAM.ann, file = paste0(out.dir, "/", Project, "-DESeq2_Resfile_collapsed_rm27_rm3_TX_SHAM_N", sigNum, ".RData"))


message("+----------                      Create results object with normCounts for collapsed samples              ------------+")

normCounts.rm <- counts(dds.collapsed.rm, normalized=TRUE)
normCounts.rm[which(rownames(normCounts.rm)=="ENSOARG00000002407"|rownames(normCounts.rm)=="ENSOARG00000012510"),]
## LEP and UCP1 counts checking

message("+ -----------------                   Run transformations                    -----------------------------------------+")

rld.collapsed.rm <- rlogTransformation(dds.collapsed.rm, blind=F)
vsd.collapsed.rm <- varianceStabilizingTransformation(dds.collapsed.rm, blind=F)

save(dds.collapsed.rm, rld.collapsed.rm, vsd.collapsed.rm, normCounts.rm,
     file = paste0(out.dir, "/", Project, "-DESeq2-dds_rld_vsd_normCounts_collapsed_rm27_rm3.RData"))

message("+--------                 DESeq2 different combinations models analysis                     --------------------------+")

ddsHTSeq.comb <-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable.new[-c(15:16,19:20),], directory=FTSeq.dir, design=~0+Comcondition)
dds_coll.comb <- collapseReplicates(ddsHTSeq.comb, groupby= ddsHTSeq.comb$Labels,  renameCols = TRUE)
dds_coll.comb <- DESeq(dds_coll.comb)
dds_coll.comb <- estimateSizeFactors(dds_coll.comb)

dds_collapsed.comb <- dds_coll.comb
normCounts.comb <- counts(dds_collapsed.comb, normalized=TRUE)

res1 <- results(dds_coll.comb, contrast=c("Comcondition", "SHAM_143", "SHAM_129"))
res2 <- results(dds_coll.comb, contrast=c("Comcondition", "TX_143", "TX_129"))
res3 <- results(dds_coll.comb, contrast=c("Comcondition", "TX_129", "SHAM_129"))
res4 <- results(dds_coll.comb, contrast=c("Comcondition", "TX_143", "SHAM_143"))

res1.df <- as.data.frame(res1)
res1.df$ensembl_gene_id <- rownames(res1.df)
res1.ann <- merge(res1.df, ensEMBL2id, by = "ensembl_gene_id")
res1.ann$description <- gsub("..Source.*", "", res1.ann$description)
res1.ann.sig <- subset(res1.ann, padj < significance & abs(log2FoldChange) >= l2fc)
print(dim(res1.ann.sig))
## 609

res2.df <- as.data.frame(res2)
res2.df$ensembl_gene_id <- rownames(res2.df)
res2.ann <- merge(res2.df, ensEMBL2id, by = "ensembl_gene_id")
res2.ann$description <- gsub("..Source.*", "", res2.ann$description)
res2.ann.sig <- subset(res2.ann, padj < significance & abs(log2FoldChange) >= l2fc)
print(dim(res2.ann.sig))
## 174


res3.df <- as.data.frame(res3)
res3.df$ensembl_gene_id <- rownames(res3.df)
res3.ann <- merge(res3.df, ensEMBL2id, by = "ensembl_gene_id")
res3.ann$description <- gsub("..Source.*", "", res3.ann$description)
res3.ann.sig <- subset(res3.ann, padj < significance & abs(log2FoldChange) >= l2fc)
print(dim(res3.ann.sig))
## 1090

res4.df <- as.data.frame(res4)
res4.df$ensembl_gene_id <- rownames(res4.df)
res4.ann <- merge(res4.df, ensEMBL2id, by = "ensembl_gene_id")
res4.ann$description <- gsub("..Source.*", "", res4.ann$description)
res4.ann.sig <- subset(res4.ann, padj < significance & abs(log2FoldChange) >= l2fc)
print(dim(res4.ann.sig))
## 1576

res5.ann <- res_TX_SHAM.ann ## compare TX vs SHAM
res5.ann.sig <- subset(res5.ann, padj < significance & abs(log2FoldChange) >= l2fc)
print(dim(res5.ann.sig))
## 1472

## add compare age only, Gestational age 143 vs 129
ddsHTSeq.age<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable.new[-c(15:16,19:20),], directory=FTSeq.dir, design=~0+GesAge)
dds_coll.age <- collapseReplicates(ddsHTSeq.age, groupby= ddsHTSeq.age$Labels,  renameCols = TRUE)
dds_coll.age <- DESeq(dds_coll.age)
dds_coll.age <- estimateSizeFactors(dds_coll.age)
normCounts.age <- counts(dds_coll.age, normalized=TRUE)

res6 <- results(dds_coll.age, contrast=c("GesAge", "143", "129"))
res6.df <- as.data.frame(res6)
res6.df$ensembl_gene_id <- rownames(res6.df)
res6.ann <- merge(res6.df, ensEMBL2id, by = "ensembl_gene_id")
res6.ann$description <- gsub("..Source.*", "", res6.ann$description)
res6.ann.sig <- subset(res6.ann, padj < significance & abs(log2FoldChange) >= l2fc)
print(dim(res6.ann.sig))
## 409

## Remove the missing log2Foldchange genes. From 27054 to 21464 genes
res1.ann <- subset(res1.ann, !is.na(log2FoldChange)); res2.ann <- subset(res2.ann, !is.na(log2FoldChange));
res3.ann <- subset(res3.ann, !is.na(log2FoldChange)); res4.ann <- subset(res4.ann, !is.na(log2FoldChange));
res5.ann <- subset(res5.ann, !is.na(log2FoldChange)); res6.ann <- subset(res6.ann, !is.na(log2FoldChange));
save(res1.ann, res2.ann, res3.ann, res4.ann, res5.ann, res6.ann, file = paste0(out.dir, "/", Project, "-SixModelCompare_Res1-6_Data.RData"))

message("+------------------ SFig 1: PCA plot for the data ------------------------------+")

sampleTBL <- sampleTable.collapsed.rm
RLD <- assay(vsd.collapsed.rm)
TOPNUM <- 500
model <- "vsd"
Group2 <- customPCA_2_shape(sampleTBL, RLD, TOPNUM, model)

message("+------------------------------------------------------------------------------+")
message("+        Check UCP1 and other individual genes as interest                     +")
message("+------------------------------------------------------------------------------+")


DDS <- dds.collapsed.rm;
CONDITIONs <- c("condition", "CombineCondition")
TITLEs1 <- paste0(c("UCP1", "DIO1", "IGF2", "LPL", "UCP2", "LEP", "IGF2R", "AKT1S1"), "_", CONDITIONs[1])
TITLEs2 <- paste0(c("UCP1", "DIO1", "IGF2", "LPL", "UCP2", "LEP", "IGF2R", "AKT1S1"), "_", CONDITIONs[2])
genes2plot <- c("ENSOARG00000012510", "ENSOARG00000006711", "ENSOARG00000003586",
               "ENSOARG00000010719", "ENSOARG00000008561", "ENSOARG00000002407",
               "ENSOARG00000004632", "ENSOARG00000010394", "ENSOARG00000013656")

for(i in 1:length(genes2plot)){
  makeGeneCountPlot(DDS, ensEMBL2id, CONDITIONs[1], TITLEs1[i], genes2plot[i], out.dir)
  makeGeneCountPlot(DDS, ensEMBL2id, CONDITIONs[2], TITLEs2[i], genes2plot[i], out.dir)
}

message("+--------------           MA plot and highlight the brown and white adipose key marker, UCP1 and LEP                   ------------+")

label_l2fc <- 3
favoritegenes <- NULL
favoritegenes <- append(favoritegenes, "ENSOARG00000012510") # UCP1
favoritegenes <- append(favoritegenes, "ENSOARG00000002407") # LEP
res.ann <- res_TX_SHAM.ann
significance <- 0.05
l2fc <- 1
res.MA <- functionCustomMAPlot(res.ann, Project, "TX_SHAM",
                               "res_TX_SHAM", significance, l2fc,
                               label_l2fc, favoritegenes, out.dir)


message("+-----------SFig 2: Volcano plot for the data, compare condition, and then split into interest groups------------------------------+")
message("+------        1) TX vs SHAM, 2) Age 129, TX vs SHAM, Age 143 TX vs SHAM; 3) TX: 143 vs 129; SHAM: 143 vs 129-           ----------+")

## load the rld file, and produce the volcano plot, TX vs SHAM
load(file = paste0(out.dir, "/", Project, "-SixModelCompare_Res1-6_Data.RData"))
resfiles.all <- list(res5.ann, res6.ann)
titles.all <- c("TX vs SHAM::N=1472", "143 vs 129::N=409")
xlabels.all <- c("log2FC (TX/SHAM)", "log2FC (143/129)")
ylabels.all <- bquote("-log"[10]~"(adj.p.value)")
Groupall.volPlot <- list()
Dnums <- c(704, 229); Unums <- c(768, 180)
for(i in 1:2){
  Groupall.volPlot[[i]] <- functionPlotDEVolcano(resfiles.all[[i]],
                                                 sig_cut=0.05,
                                                 logfc_cut=1,
                                                 title=titles.all[i],
                                                 xrange=c(-7.5,7.5,2.5),
                                                 yrange=c(0,80,20),
                                                 topN=10,
                                                 xlabel=xlabels.all[i],
                                                 ylabel=ylabels.all,
                                                 dnum=Dnums[i],
                                                 unum=Unums[i])
}

pdf(paste0(out.dir, "/", Project, "-ModelComparison_Group56_VolcanoPlot_23_07_2019.pdf"), height = 8, width = 16)
plot_grid(Groupall.volPlot[[1]], Groupall.volPlot[[2]], ncol=2)
dev.off()

## defined groups within treatment or within gestational age.
titles <- c( "SHAM (143 vs 129)::N=609", "TX (143 vs 129)::N=174","129 (TX vs SHAM)::N=1090", "143 (TX vs SHAM)::N=1576")
resfiles <- list(res1.ann, res2.ann,res3.ann, res4.ann)
xlabels <- rep(xlabels.all, each = 2)
ylabels <- ylabels.all
Group.volPlot <- list()
Dnums <- c(377,88,465,707); Unums <- c(232,86,625,869)

for(i in 1:4){
  Group.volPlot[[i]] <- functionPlotDEVolcano(resfiles[[i]],
                                                 sig_cut=0.05,
                                                 logfc_cut=1,
                                                 title=titles[i],
                                                 xrange=c(-7.5,7.5,2.5),
                                                 yrange=c(0,80,20),
                                                 topN=10,
                                                 xlabel=xlabels[i],
                                                 ylabel=ylabels,
                                                 dnum=Dnums[i],
                                                 unum=Unums[i])
}

pdf(paste0(out.dir, "/", Project, "-ModelComparison_Group4_VolcanoPlot_23_07_2019.pdf"), height = 10, width = 10)
plot_grid(Group.volPlot[[1]], Group.volPlot[[2]], Group.volPlot[[3]], Group.volPlot[[4]], nrow=2)
dev.off()

message("+-----        Heatmap for the top 262 genes with cutoff abs(l2fc) >= 2 and padj < 0.05      -------------------+")

houtfile1 <- paste0(out.dir, "/", Project, "-Heatmap_sortData_padj0.05_l2fc2_top272_rm27_rm3_23_07_2019.RData")
select.no <- 272
datasort <- heatmap_data_fn(res_TX_SHAM.ann, psig.level=0.05, lfc.level=2, select.no, houtfile1)

load(houtfile1);
s12.sig <- res.flc.padj;
sampleTable.collapsed.rm$sample_id <- paste0(sampleTable.collapsed.rm$condition, "_", sampleTable.collapsed.rm$GestationalAge,
                                          "_", sampleTable.collapsed.rm$Labels)
# check gestational age as a factor
outHfile1 <- paste0(out.dir, "/", Project, "-Heatmap_TX_SHAM_padj0.05_lfc2_top272_rm27_rm3_rmLC_23_07_2019.pdf")
rlds <- vsd.collapsed.rm
genes2plot <- s12.sig
ensE <- ensEMBL2id
samTable <- sampleTable.collapsed.rm
outfile <- outHfile1
width <- 7; height = 10
heatmap_mat_fn_age(rlds, genes2plot, ensE, samTable, outfile, width, height)

# check the number of genes list with padj < 0.05 & abs(l2fc) >=3
nums <- dim(subset(res_TX_SHAM.ann, padj<0.05&abs(log2FoldChange)>=3))[1]
houtfile2 <- paste0(out.dir, "/", Project, "-Heatmap_sortData_padj0.05_l2fc3_top", nums, "_rm27_rm3_rmLC_23_07_2019.RData")
select.no <- nums
datasort <- heatmap_data_fn(res_TX_SHAM.ann, psig.level=0.05, lfc.level=3, select.no, houtfile2)

load(houtfile2);
s12.sig <- res.flc.padj;

outHfile2 <- paste0(out.dir, "/", Project, "-Heatmap_TX_SHAM_padj0.05_lfc3_top", nums,"_rm27_rm3_rmLC_23_07_2019.pdf")
heatmap_mat_fn_age(rlds, s12.sig, ensEMBL2id, samTable, outHfile2, width=7, height=10)

message("+----------     Generate differential compare differetial analysis, 143 TXvsSHAM compare to 129 TXvsSHAM     ----------+")
message("+----------     Generate differential compare differetial analysis, TX 143vs129 compare to SHAM 143vs129     ----------+")

clustA <- res4.ann; clustB <- res3.ann;
significance1 <- 0.05
logfc_cut1 <- 1
topN <- 10
significance2 <- 0.01
logfc_cut2 <- 1.5
nameTB <- "Table.DEvsDE.143_TX_vs_143_SHAM.VS.129_TX_vs_129_SHAM_blue_N"
nameTG <- "Table.DEvsDE.143_TX_vs_143_SHAM.VS.129_TX_vs_129_SHAM_green_N"
nameTP <- "Table.DEvsDE.143_TX_vs_143_SHAM.VS.129_TX_vs_129_SHAM_purple_N"
nameTY <- "Table.DEvsDE.143_TX_vs_143_SHAM.VS.129_TX_vs_129_SHAM_grey_N"

clustA.lab <- "143 (TX/SHAM)"; clustB.lab <- "129 (TX/SHAM)"
pdf(paste0(out.dir, "/", Project, "-CorrPlot_143_TXvSHAM_vs_129_TXvSHAM_p0.05_l2fc1_p0.01_l2fc1.5_31_07_2019.pdf"), height =  10, width = 10)
corrD.plot1 <- functionPlotDECorrelation(clustA,clustB,  clustA.lab, clustB.lab, significance1, logfc_cut1, topN, out.dir, nameTB, nameTG, nameTP, nameTY)
print(corrD.plot1)
corrD.plot2 <- functionPlotDECorrelation(clustA,clustB,  clustA.lab, clustB.lab, significance2, logfc_cut2, topN, out.dir, nameTB, nameTG, nameTP, nameTY)
print(corrD.plot2)
dev.off()

clustA <- res2.ann; clustB <- res1.ann;
nameTB <- "Table.DEvsDE.TX_143_vs_129.VS.SHAM_143_vs_129_blue_N"
nameTG <- "Table.DEvsDE.TX_143_vs_129.VS.SHAM_143_vs_129_green_N"
nameTP <- "Table.DEvsDE.TX_143_vs_129.VS.SHAM_143_vs_129_purple_N"
nameTY <- "Table.DEvsDE.TX_143_vs_129.VS.SHAM_143_vs_129_grey_N"

clustA.lab <- "TX (143/129)"; clustB.lab <- "SHAM (143/129)"
pdf(paste0(out.dir, "/", Project, "-CorrPlot_TX_143v129_vs_SHAM_143v129_p0.05_l2fc1_p0.01_l2fc1.5_31_07_2019.pdf"), height =  10, width = 10)
corrD.plot1 <- functionPlotDECorrelation(clustA,clustB, clustA.lab, clustB.lab, significance1, logfc_cut1, topN, out.dir, nameTB, nameTG, nameTP, nameTY)
print(corrD.plot1)
corrD.plot2 <- functionPlotDECorrelation(clustA,clustB,  clustA.lab, clustB.lab, significance2, logfc_cut2, topN, out.dir, nameTB, nameTG, nameTP, nameTY)
print(corrD.plot2)
dev.off()

message("+-----                          Key markers for brown/beige and white adipocytes to check         ------------+")

## paper Nat Rev Mol Cell Biol. 2016 November ; 17(11): 691–702. doi:10.1038/nrm.2016.96,
## identified mouse adipose key markers as below, In Table 1.

BrownORBeige.markers <- c("UCP1", "DIO2", "CIDEA", "PPARGC1A", "PPARA", "COX7A1", "PRDM16", "EBF2")
## EBF2 is both adipocyte and Pred.
## COX8B is excluded as in sheep which is unannoted.

Brown.markers <- c("LHX8", "PDK4")
## "ZIC1" has no mapping in our data, but exist in gtf file.
## "EVA1" and  "EPSTL1" no annotations.

Beige.markers <- c("CITED1", "SHOX2", "TMEM26")
## "TBX1", no mapping reads.
## "CD137", "PAT2", "P2RX5" no annotations.

White.markers <- c("LEP", "RETN", "AGT", "WT1")
## WT1 is pred.

BrownBeigeWhite <- c("ADIPOQ", "FABP4",  "PDGFRA", "CD34", "PPARG")
## CD34 and PDGFRA are pred, give different line colors.
## "CEBPB", "SCA1", "PREF1", "CD29" no annotations

# BrownTransition <- read.xlsx(paste0(Base.dir, "/Paper_Basse_BMC_Genomics_2015/12864_2015_1405_MOESM5_ESM.xlsx"), sheetIndex=1)
# TransitionWhite <- read.xlsx(paste0(Base.dir, "/Paper_Basse_BMC_Genomics_2015/12864_2015_1405_MOESM6_ESM.xlsx"), sheetIndex=1)

## Use the normCounts.rm matrix to produce the line/dot plot across four group, SHAM_129, SHAM_143, TX_129, TX_143
normCounts.rm.dat <- as.data.frame(normCounts.rm)
normCounts.rm.dat$ensembl_gene_id <- rownames(normCounts.rm.dat)
normCounts.rm.dat.annot <- merge(normCounts.rm.dat, ensEMBL2id, by = "ensembl_gene_id")

allmarkers <- cbind(c(Brown.markers, Beige.markers, White.markers, BrownORBeige.markers, BrownBeigeWhite),
                    c(rep(c(1:5), c(2,3,4,8,5))))
colnames(allmarkers) <- c("external_gene_name", "group")

select.subNormCounts <- normCounts.rm.dat.annot[normCounts.rm.dat.annot[,37]%in%allmarkers[,1],]
colnames(select.subNormCounts) <- c("ensembl_gene_id", sampleTable.collapsed.rm27.rm3$Comcondition, "external_gene_name",
                                    "entezgene_id", "description")
SHAM_129.cols <- which(grepl("SHAM_129", colnames(select.subNormCounts))==T)
TX_129.cols <- which(grepl("TX_129", colnames(select.subNormCounts))==T)
SHAM_143.cols <- which(grepl("SHAM_143", colnames(select.subNormCounts))==T)
TX_143.cols <- which(grepl("TX_143", colnames(select.subNormCounts))==T)

meanSHAM129 <- rowMeans(select.subNormCounts[,SHAM_129.cols]);
meanSHAM143 <- rowMeans(select.subNormCounts[,SHAM_143.cols]);
meanTX129 <- rowMeans(select.subNormCounts[,TX_129.cols]);
meanTX143 <- rowMeans(select.subNormCounts[,TX_143.cols]);
IDs <- rep(c("Brown/Beige", "White", "White", "Beige", "Brown/Beige", "White", "Brown", "Beige", "Brown/Beige", "Brown/Beige",
             "BrownBeigeWhite", "Brown/Beige", "Brown/Beige", "Brown", "BrownBeigeWhite", "White", "BrownBeigeWhite",
             "Brown/Beige", "Beige", "BrownBeigeWhite", "Brown/Beige", "BrownBeigeWhite"), each=4)
genes <- as.character(rep(select.subNormCounts$external_gene_name, each=4))
groups <- rep(c("SHAM_129", "SHAM_143", "TX_129", "TX_143"), length=88)

Counts <- NULL
for( i in 1:22){
  sCounts <- c(meanSHAM129[i], meanSHAM143[i], meanTX129[i], meanTX143[i])/meanSHAM129[i]
  Counts <- c(Counts, sCounts)
  Counts
}

newPlotdat <- tibble(genes, groups, IDs, Counts)

pdf(paste0(out.dir, "/", Project, "-BubblePlot_select_Brown_Beige_White_Mouse_01_08_2019.pdf"))
ggplot(newPlotdat, aes(y =genes, x = groups), alpha=0.75) +
  geom_point(aes(size = Counts, colour = IDs)) +
  guides(size=guide_legend(title="Normalised Counts Ratio to SHAM 129")) +
  guides(colour = guide_legend(title = "Adipose Key Markers Category")) +
  xlab("") + ylab("") +
  theme_bw()

dev.off()



message("+--------------Produce a bar plot for 6 model comparisons up/down regulated genes, use padj < 0.05 & l2fc >=1 -----+")

go.significance <- 0.05
go.l2fc <- 1
geneL_TX_SHAM <- subset(res5.ann, abs(log2FoldChange)>=go.l2fc & padj < go.significance);
geneL_TX_143vs129 <- subset(res2.ann, abs(log2FoldChange)>=go.l2fc & padj < go.significance)
geneL_SHAM_143vs129 <- subset(res1.ann, abs(log2FoldChange)>=go.l2fc & padj < go.significance)
geneL_129_TXvsSHAM <- subset(res3.ann, abs(log2FoldChange)>=go.l2fc & padj < go.significance)
geneL_143_TXvsSHAM <- subset(res4.ann, abs(log2FoldChange)>=go.l2fc & padj < go.significance)
geneL_143vs129 <- subset(res6.ann, abs(log2FoldChange)>=go.l2fc & padj < go.significance)

Models <- rep(c("Ges 143: TX vs SHAM", "TX vs SHAM",  "Ges 129: TX vs SHAM", "SHAM: 143 vs 129", "Ges 143 vs Ges 129", "TX: 143 vs 129"),
              length = 12)
DEGs <- c(-sum(geneL_143_TXvsSHAM$log2FoldChange<0), -sum(geneL_TX_SHAM$log2FoldChange<0),
          -sum(geneL_129_TXvsSHAM$log2FoldChange<0), -sum(geneL_SHAM_143vs129$log2FoldChange<0),
          -sum(geneL_143vs129$log2FoldChange<0), -sum(geneL_TX_143vs129$log2FoldChange<0),
          sum(geneL_143_TXvsSHAM$log2FoldChange >0), sum(geneL_TX_SHAM$log2FoldChange >0),
          sum(geneL_129_TXvsSHAM$log2FoldChange >0), sum(geneL_SHAM_143vs129$log2FoldChange >0),
          sum(geneL_143vs129$log2FoldChange >0), sum(geneL_TX_143vs129$log2FoldChange >0) )
Nums <- c(sum(geneL_143_TXvsSHAM$log2FoldChange<0), sum(geneL_TX_SHAM$log2FoldChange<0),
          sum(geneL_129_TXvsSHAM$log2FoldChange<0), sum(geneL_SHAM_143vs129$log2FoldChange<0),
          sum(geneL_143vs129$log2FoldChange<0), sum(geneL_TX_143vs129$log2FoldChange<0),
          sum(geneL_143_TXvsSHAM$log2FoldChange >0), sum(geneL_TX_SHAM$log2FoldChange >0),
          sum(geneL_129_TXvsSHAM$log2FoldChange >0), sum(geneL_SHAM_143vs129$log2FoldChange >0),
          sum(geneL_143vs129$log2FoldChange >0), sum(geneL_TX_143vs129$log2FoldChange >0) )
Totals <-  rep(c(1576, 1472, 1090, 609, 409, 174), length=12)

Direc <- rep(c("Down", "Up"), each=6)
barPlotdat <- data.frame(Models, Direc, Totals, DEGs, Nums)

brks <- seq(-900, 900, 300)
lbls = as.character(c(seq(900, 0, -300), seq(300, 900, 300)))
pdf(paste0(out.dir, "/", Project, "-Barplot_SigDEGs_summary_Model1-6_01_08_2019.pdf"))
barplotSum <- ggplot(barPlotdat, aes(x = reorder(Models, -Totals), y = DEGs, fill = Direc)) +   # Fill column
                              geom_bar(stat = "identity", width = 0.5) +   # draw the bars
                              #geom_text(aes(label = Nums), size = 4, hjust = 0.25, vjust = 0.8) +
                              scale_y_continuous(breaks = brks,   # Breaks
                                                 labels = lbls) + # Labels
                              coord_flip() +  # Flip axes
                              labs(title="Number of Up/Down Regulated Genes summary") +
                              xlab("") +
                              ylab("") +
                              theme_tufte() +  # Tufte theme from ggfortify
                              theme(plot.title = element_text(hjust = .2),
                                    axis.ticks = element_blank()) +   # Centre plot title
                               scale_fill_manual( values = c("blue", "red"))  # Color palette
print(barplotSum)
dev.off()

message("+------- 6 models upsetR plot 1) all intersection groups 2) only show the number of unique groups after intersection-+")

UpsetR.datList <- list(geneL_TX_SHAM$ensembl_gene_id, geneL_TX_143vs129$ensembl_gene_id, geneL_SHAM_143vs129$ensembl_gene_id,
                       geneL_129_TXvsSHAM$ensembl_gene_id, geneL_143_TXvsSHAM$ensembl_gene_id, geneL_143vs129$ensembl_gene_id)
names(UpsetR.datList) <- c("TX vs SHAM (1472)", "TX_143 vs TX 129 (174)", "SHAM 143 vs SHAM 129 (609)",
                           "129 TX vs 129 SHAM (1090)", "143 TX vs 143 SHAM (1576)", "143 vs 129 (409)")

pdf(paste0(out.dir, "/", Project, "-DESeq2_6Comparisons_sigGenes_UpsetR_all_05_08_2019.pdf"), height = 5, width = 10)
upset(fromList(UpsetR.datList), order.by = "freq", nsets = 6, nintersects = 47, sets.bar.color = "#56B4E9")
dev.off()


pdf(paste0(out.dir, "/", Project, "-DESeq2_6Comparisons_sigGenes_UpsetR_onlyUnique_05_08_2019.pdf"), height = 5, width = 10)
upset(fromList(UpsetR.datList), order.by = c("freq", "degree"), nsets = 6, nintersects = 6, sets.bar.color = "#56B4E9", empty.intersections = "on", keep.order = F)
dev.off()


message("+-----------------------------------------------------------------------------------------------------------------+")
message("+--------------------                         Gene Ontology analysis                              ----------------+")
message("+-----------------------------------------------------------------------------------------------------------------+")
## Basic libraries calling

suppressPackageStartupMessages({
library('clusterProfiler')
library('DOSE')
library('GSEABase')
library('GEOquery')
library('UniProt.ws')
library('rgdal')
library('AnnotationDbi')
library('org.Bt.eg.db')
library('GSEABase')
library('enrichplot')
library('pathview')
library('AnnotationHub')
library('circlize')
library(ensembldb)
})

hub <- AnnotationHub()
## snapshotDate(): 2018-10-24
query(hub, "ovis_aries")

##  AnnotationHub with 151 records
###  snapshotDate(): 2018-10-24
###  $dataprovider: Ensembl, ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
###  $species: Ovis aries
###  $rdataclass: FaFile, TwoBitFile, GRanges, OrgDb
###  additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
###    rdatapath, sourceurl, sourcetype
###  retrieve records with, e.g., 'object[["AH8773"]]'

##  title
##  AH8773  | Ovis_aries.Oar_v3.1.74.gtf
##  AH10704 | Ovis_aries.Oar_v3.1.75.gtf
##  AH10995 | Ovis_aries.Oar_v3.1.75.cdna.all.fa
##  AH10996 | Ovis_aries.Oar_v3.1.75.dna_rm.toplevel.fa
##  AH10997 | Ovis_aries.Oar_v3.1.75.dna_sm.toplevel.fa
##  ...       ...
##  AH65974 | Ovis_aries.Oar_v3.1.cdna.all.2bit
##  AH65975 | Ovis_aries.Oar_v3.1.dna_rm.toplevel.2bit
##  AH65976 | Ovis_aries.Oar_v3.1.dna_sm.toplevel.2bit
##  AH65977 | Ovis_aries.Oar_v3.1.ncrna.2bit
##  AH66408 | org.Ovis_aries.eg.sqlite

Sheep <- hub[["AH66408"]]


message("+-----------     Get a list of homolog genes between sheep and cow, most close species,     ----------------------- +")
message("+------------- Due to the annotation of sheep has lots of missing entrezid, use cow as a ref -----------------------+")

sheep     =  useEnsembl(biomart="ensembl", dataset="oaries_gene_ensembl")
cow       =  useEnsembl(biomart="ensembl", dataset="btaurus_gene_ensembl")
sheep.ens <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "entrezgene_id"),mart = sheep)

message("+-----  Apply sig genes or homolog genes Entrez ID do GO analysis with cluterProfiler           -----------------------+")
message("+----- Due to the poor annotation of sheep, we also perform the homolog Cow GO analysis         -----------------------+")

## generate the data input for GO clusterProfiler. choose the entrezid, external_gene_name and log2foldchange.

select.cols <- c(1,8,9,3)
geneL_list <- list(geneL_TX_SHAM, geneL_TX_143vs129, geneL_SHAM_143vs129, geneL_129_TXvsSHAM, geneL_143_TXvsSHAM, geneL_143vs129)
genes_GO_List <- list(); genes_KEGG_List <- list()
BP.sheep <- list(); MF.sheep <- list(); CC.sheep <- list(); KEGG.sheep <- list()
BP.sheep.ori <- list(); KEGG.sheep.ori <- list()

for(i in 1:6){
      genes_GO_List[[i]] <- geneL_list[[i]][order(-geneL_list[[i]]$log2FoldChange), select.cols]
      colnames(genes_GO_List[[i]]) <- c("ENSEMBL", "SYMBOL", "ENTREZ", "L2FC")
      gc()
      genes_KEGG <- genes_GO_List[[i]][!is.na(genes_GO_List[[i]]$ENTREZ),]
      genes_KEGGs <- genes_KEGG[,4]
      names(genes_KEGGs) <- genes_KEGG[,3]
      genes_KEGG_List[[i]] <- genes_KEGGs
      ##

      BP.sheep[[i]] <- enrichGO(genes_GO_List[[i]]$ENTREZ,
                                OrgDb = Sheep,
                                ont="BP",
                                minGSSize = 5,
                                keyType = "ENTREZID",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = FALSE)
      BP.sheep.ori[[i]] <- BP.sheep[[i]]
      BP.sheep[[i]] <- setReadable(BP.sheep[[i]], keyType = "ENTREZID", OrgDb = Sheep)
      ##
      CC.sheep[[i]] <- enrichGO(genes_GO_List[[i]]$ENTREZ,
                                OrgDb = Sheep,
                                ont="CC",
                                minGSSize = 5,
                                keyType = "ENTREZID",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
      CC.sheep[[i]] <- setReadable(CC.sheep[[i]], keyType = "ENTREZID", OrgDb = Sheep)
      ##
      MF.sheep[[i]] <- enrichGO(genes_GO_List[[i]]$ENTREZ,
                                OrgDb = Sheep,
                                ont="MF",
                                minGSSize = 5,
                                keyType = "ENTREZID",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                readable = TRUE)
      MF.sheep[[i]] <- setReadable(MF.sheep[[i]], keyType = "ENTREZID", OrgDb = Sheep)
      ##
      KEGG.sheep[[i]] <- enrichKEGG(names(genes_KEGG_List[[i]]),
                                organism = 'oas',
                                minGSSize = 5,
                                keyType = "kegg",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05)
      KEGG.sheep.ori[[i]] <- KEGG.sheep[[i]]
      KEGG.sheep[[i]] <- setReadable(KEGG.sheep[[i]], keyType = "ENTREZID", OrgDb = Sheep)

}

save(BP.sheep.ori, KEGG.sheep.ori, BP.sheep, KEGG.sheep, file = paste0(out.dir, "/", Project, "-GO_KEGG_pathway_p0.05_l2fc1_q0.05_08_08_2019.RData"))

message("+-------------Produce the GO/KEGG heatmap plot for six models common/unique pathways ------------+")

## KEGG pathways : model TX vs SHAM (1140): 44 pathways; model TX 143 vs 129 (147): 1 pathway; model SHAM 143 vs 129 (504): 2 pathways;
##                 model 129 TX vs SHAM (882): 47 pathways; model 143 TX vs SHAM (1230): 29 pathways; model 143 vs 129 (327): 1 pathway.

KEGG.pathways <- unique(unlist(lapply(KEGG.sheep, function(x) x$Description)))
KEGG.pathways.datList <- lapply(KEGG.sheep, function(x) data.frame(x$Description, x$p.adjust))
KEGG.pathway.dat <- Reduce(function(x,y) merge(x,y,by="x.Description",all=TRUE),KEGG.pathways.datList)
KEGG.pathway.dat <- KEGG.pathway.dat[,c(1,2,5,6,4,3,7)]
colnames(KEGG.pathway.dat) <- c("Pathway", "TXvsSHAM", "129:TXvsSHAM", "143:TXvsSHAM", "SHAM:143vs129",  "TX:143vs129", "143vs129")
KEGG.pathway.final.dat <- KEGG.pathway.dat[,-1]
rownames(KEGG.pathway.final.dat) <- KEGG.pathway.dat[,1]
#KEGG.pathway.final.dat[is.na(KEGG.pathway.final.dat)] <- 0
KEGG.pathway.final.mat <- as.matrix(KEGG.pathway.final.dat)

## use ggplot2 to produce heatmap
KEGG.gplotdat <- KEGG.pathway.dat
#KEGG.gplotdat[is.na(KEGG.gplotdat)] <- 0
KEGG.gplotdat$Pathway <- with(KEGG.gplotdat, reorder(Pathway, TXvsSHAM))
KEGG.gplotm <- melt(KEGG.gplotdat)
KEGG.gplotm <- ddply(KEGG.gplotm, .(variable), transform)

pdf(paste0(out.dir, "/", Project, "-KEGGpathway_6models_Heatmap_05_08_2019.pdf"), height=10, width = 5)
base_size <- 9
ggplot(KEGG.gplotm, aes(variable, Pathway)) +
      geom_tile(aes(fill = value), colour = "grey") +
      scale_fill_gradient(low = "red", high = "steelblue3", guide="colorbar",na.value="white", name= "p.adj") +
      theme_grey(base_size = base_size) +
      labs(x = "", y = "") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(fill = "values")
dev.off()



## use ComplexHeatmap to produce the heatmap.

col_fun = colorRamp2(c(0, 0.01, 0.05), c("red", "yellow", "green"))
#col_fun(seq(-3, 3))

KEGG.ht <- Heatmap(KEGG.pathway.final.mat,
        na_col = "grey",
        row_names_side = "left",
        cluster_rows=F,
        cluster_row_slices = F,
        cluster_columns = F,
        cluster_column_slices = F,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 9),
        col = col_fun,
        width = 16,
        height = 8,
        heatmap_legend_param = list(direction = "vertical", title="p.adj"),

        column_order = sort(colnames(KEGG.pathway.final.mat))
)
draw(KEGG.ht, heatmap_legend_side="right")

message("+---------------                     Biological Process for six models             ----------------------------+")

## Biological Process : model TX vs SHAM (1140): 24 pathways; model TX 143 vs 129 (147): 0 pathway; model SHAM 143 vs 129 (504): 2 pathways;
##                 model 129 TX vs SHAM (882): 100 pathways; model 143 TX vs SHAM (1230): 7 pathways; model 143 vs 129 (327): 0 pathway.

BP.pathways <- unique(unlist(lapply(BP.sheep, function(x) x$Description)))
BP.pathways.datList <- lapply(BP.sheep, function(x) data.frame(x$Description, x$p.adjust))
BP.pathway.dat <- Reduce(function(x,y) merge(x,y,by="x.Description",all=TRUE),BP.pathways.datList)
BP.pathway.dat <- BP.pathway.dat[,c(1,2,5,6,4,3,7)]
colnames(BP.pathway.dat) <- c("Pathway", "TXvsSHAM", "129:TXvsSHAM", "143:TXvsSHAM", "SHAM:143vs129",  "TX:143vs129", "143vs129")

## use ggplot2 to produce heatmap
BP.gplotdat <- BP.pathway.dat
BP.gplotdat$Pathway <- with(BP.gplotdat, reorder(Pathway, TXvsSHAM))
BP.gplotm <- melt(BP.gplotdat)
BP.gplotm <- ddply(BP.gplotm, .(variable), transform)

pdf(paste0(out.dir, "/", Project, "-BPpathway_4models_Heatmap_05_08_2019.pdf"), height=12, width = 5)
base_size <- 9
ggplot(BP.gplotm, aes(variable, Pathway)) +
      geom_tile(aes(fill = value), colour = "grey") +
      scale_fill_gradient(low = "red", high = "steelblue3", guide="colorbar",na.value="white", name= "p.adj") +
      theme_grey(base_size = base_size) +
      labs(x = "", y = "") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(fill = "values")
dev.off()


message("+-----       Produce the barplot with the number of up/down regulated genes and zscore, and padj  ----------------+")

KEGG_list_up <- list()
KEGG.sheep.TS <- as.data.frame(KEGG.sheep.ori[[1]])

for (i in 1:nrow(KEGG.sheep.TS)){
  df_tmp <- res5.ann[res5.ann$entrezgene_id %in% unlist(strsplit(KEGG.sheep.TS[i, "geneID"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  KEGG_list_up[[i]] <-tmp_up
}

KEGG.sheep.TS$UP <- as.numeric(KEGG_list_up)
KEGG.sheep.TS$DOWN <- KEGG.sheep.TS$Count - KEGG.sheep.TS$UP
KEGG.sheep.TS$zscore <- (KEGG.sheep.TS$UP-KEGG.sheep.TS$DOWN)/sqrt(KEGG.sheep.TS$Count)
KEGG.sheep.TS$DOWN <- -KEGG.sheep.TS$DOWN


KEGG.sheep.TS_molten <- melt(KEGG.sheep.TS[,c(2,6:7,9:11)],
                      id.vars=c("Description","p.adjust","qvalue","Count") )

brks <- seq(-30, 30, 10)
lbls = as.character(c(seq(30, 0, -10), seq(10, 30, 10)))
zscore <- round(KEGG.sheep.TS$zscore, digits = 3)

pdf(paste0(out.dir, "/", Project, "-Kegg_sheep_TXvsSHAM_GeneCount_Barplot_N44_08_08_19.pdf"), width=12, height=8)
barplotSum <- ggplot(KEGG.sheep.TS_molten, aes(x = reorder(Description, -qvalue), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = 0.8) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls,
                     limits=c(-30,30)) + # Labels

  coord_flip() +  # Flip axes
  labs(title="KEGG Pathway") +
  xlab("") +
  ylab("") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = .2),
        axis.ticks = element_blank()) +   # Centre plot title
  scale_fill_manual(values = c("red", "blue"))+  # Color palette
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  annotate("text", x = 1 , y =25, label = zscore[44], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 2 , y =25, label = zscore[43], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 3 , y =25, label = zscore[42], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 4 , y =25, label = zscore[41], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 5 , y =25, label = zscore[40], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 6 , y =25, label = zscore[39], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 7 , y =25, label = zscore[38], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 8 , y =25, label = zscore[37], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 9 , y =25, label = zscore[36], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 10 , y =25, label = zscore[35], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 11 , y =25, label = zscore[34], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 12 , y =25, label = zscore[33], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 13 , y =25, label = zscore[32], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 14 , y =25, label = zscore[31], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 15 , y =25, label = zscore[30], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 16 , y =25, label = zscore[29], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 17 , y =25, label = zscore[28], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 18 , y =25, label = zscore[27], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 19 , y =25, label = zscore[26], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 20 , y =25, label = zscore[25], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 21 , y =25, label = zscore[24], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 22 , y =25, label = zscore[23], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 23 , y =25, label = zscore[22], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 24 , y =25, label = zscore[21], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 25 , y =25, label = zscore[20], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 26 , y =25, label = zscore[19], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 27 , y =25, label = zscore[18], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 28 , y =25, label = zscore[17], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 29 , y =25, label = zscore[16], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 30 , y =25, label = zscore[15], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 31 , y =25, label = zscore[14], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 32 , y =25, label = zscore[13], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 33 , y =25, label = zscore[12], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 34 , y =25, label = zscore[11], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 35 , y =25, label = zscore[10], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 36 , y =25, label = zscore[9], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 37 , y =25, label = zscore[8], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 38 , y =25, label = zscore[7], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 39 , y =25, label = zscore[6], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 40 , y =25, label = zscore[5], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 41 , y =25, label = zscore[4], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 42 , y =25, label = zscore[3], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 43 , y =25, label = zscore[2], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 44 , y =25, label = zscore[1], vjust = 0.5, hjust = 0.5)


barplotSum

dev.off()

## all the biological process bar plot for TX vs SHAM ##

BP_list_up <- list()
BP.sheep.TS <- as.data.frame(BP.sheep.ori[[1]])

for (i in 1:nrow(BP.sheep.TS)){
  df_tmp <- res5.ann[res5.ann$entrezgene_id %in% unlist(strsplit(BP.sheep.TS[i, "geneID"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  BP_list_up[[i]] <-tmp_up
}

BP.sheep.TS$UP <- as.numeric(BP_list_up)
BP.sheep.TS$DOWN <- BP.sheep.TS$Count - BP.sheep.TS$UP
BP.sheep.TS$zscore <- (BP.sheep.TS$UP-BP.sheep.TS$DOWN)/sqrt(BP.sheep.TS$Count)
BP.sheep.TS$DOWN <- -BP.sheep.TS$DOWN


BP.sheep.TS_molten <- melt(BP.sheep.TS[,c(2,6:7,9:11)],
                      id.vars=c("Description","p.adjust","qvalue","Count") )

brks <- seq(-40, 40, 10)
lbls = as.character(c(seq(40, 0, -10), seq(10, 40, 10)))
zscore <- round(BP.sheep.TS$zscore, digits = 3)

pdf(paste0(out.dir, "/", Project, "-BP_sheep_TXvsSHAM_GeneCount_Barplot_N24_08_08_19.pdf"), width=12, height=8)
barplotSum1 <- ggplot(BP.sheep.TS_molten, aes(x = reorder(Description, -qvalue), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = 0.8) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls,
                     limits=c(-40,40)) + # Labels

  coord_flip() +  # Flip axes
  labs(title="Biological Process") +
  xlab("") +
  ylab("") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = .2),
        axis.ticks = element_blank()) +   # Centre plot title
  scale_fill_manual(values = c("red", "blue"))+  # Color palette
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  annotate("text", x = 1 , y =35, label = zscore[24], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 2 , y =35, label = zscore[23], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 3 , y =35, label = zscore[22], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 4 , y =35, label = zscore[21], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 5 , y =35, label = zscore[20], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 6 , y =35, label = zscore[19], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 7 , y =35, label = zscore[18], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 8 , y =35, label = zscore[17], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 9 , y =35, label = zscore[16], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 10 , y =35, label = zscore[15], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 11 , y =35, label = zscore[14], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 12 , y =35, label = zscore[13], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 13 , y =35, label = zscore[12], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 14 , y =35, label = zscore[11], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 15 , y =35, label = zscore[10], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 16 , y =35, label = zscore[9], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 17 , y =35, label = zscore[8], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 18 , y =35, label = zscore[7], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 19 , y =35, label = zscore[6], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 20 , y =35, label = zscore[5], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 21 , y =35, label = zscore[4], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 22 , y =35, label = zscore[3], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 23 , y =35, label = zscore[2], vjust = 0.5, hjust = 0.5) +
  annotate("text", x = 24 , y =35, label = zscore[1], vjust = 0.5, hjust = 0.5)


barplotSum1

dev.off()
