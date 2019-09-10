## !/usr/local/bin/R
## Rscript for the functions that will call for the RNASeq analysis of sheep Adipose tissue
## Author: Xiaohui Zhao (xz289@cam.ac.uk)
## revised on 23/07/2019

message("+----- Customised PCA plot function, colored by TX and SHAM, shaped by treatment+gestational age -----+")
customPCA_2_shape <- function(sampleTBL, RLD, TOPNUM, model) {

  rv     <- rowVars(RLD)
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD[select, ]))

  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

  scores    <- data.frame(sampleName=sampleTBL$sampleName, pca$x, condition=sampleTBL$condition,
                          GesAge=sampleTBL$GestationalAge, Comcondition=sampleTBL$CombineCondition)

  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, shape=Comcondition,
                                label=scores$sampleName) ) +

    geom_point(size = 3, alpha=0.75 ) +
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    xlab(pc1lab) + ylab(pc2lab) + #coord_fixed() +
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    scale_shape_manual(name="Comcondition", values = c(1,2, 19, 24)) +
    theme(text = element_text(size=elementTextSize))

  pdf(paste0(out.dir, "/", Project, "-ColorPCA_shapeAge_TX_SHAM_rm27_rm3_collapsed_newshape_", model, "_23_07_19.pdf"),
  width= 8 , height = 6)
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition,
                                   shape=Comcondition, label=scores$sampleName )) +
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(colour=condition, fill=condition, group = condition)) +
    geom_point(size = 3, alpha=0.75 ) +
    xlab(pc1lab) + ylab(pc2lab) +
    scale_shape_manual(name="CombineCondition", values = c(1,2, 19, 17)) +
    scale_color_manual(name="Treatment", values = c("SHAM"="red", "TX"= "blue"))  +
    scale_fill_manual(name="Treatment", values=c("SHAM"="red", "TX"= "blue")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize))
  print(plt.pca.nl)
  dev.off()

  loadings                 <- as.data.frame(pca$rotation)
  loadings$ensembl_gene_id <- rownames(loadings)
  loadings                 <- merge(loadings, ensEMBL2id, by="ensembl_gene_id")

  pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
  pca.1.name    <- subset(pca.1, external_gene_name!="")
  pca.1.25      <- pca.1.name[c(1:25),]

  pdf(paste0(out.dir, "/", Project, "-ColorPCA_PC1_TX_SHAM_rm27_rm3_collapsed_newshape_", model, "_23_07_19.pdf"), height=5, width=8)
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(pca.1.25$external_gene_name,levels=unique(pca.1.25$external_gene_name)), y=PC1)) +
    geom_point(size = 3 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  print(pca.1.25.plot)
  dev.off()

  pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
  pca.2.name    <- subset(pca.2, external_gene_name!="")
  pca.2.25      <- pca.2.name[c(1:25),]
  pdf(paste0(out.dir, "/", Project, "-ColorPCA_PC2_TX_SHAM_rm27_rm3_collapsed_newshape_", model, "_23_07_19.pdf"), height = 5, width = 8)
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(pca.2.25$external_gene_name,levels=unique(pca.2.25$external_gene_name)), y=PC2)) +
    geom_point(size = 3 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  print(pca.2.25.plot)
  dev.off()
}
message("+------     Individual genes counts plot as interested by different conditions settings --------------+")
makeGeneCountPlot <- function(DDS, ensEMBL2id, CONDITION, TITLE, gene2plot,outdir) {
  #
  # Plot the normalised read counts for a specified gene
  #
  if(missing(outdir)){ outdir = "" }
  else( outdir <- paste(outdir, "/", sep=""))

  genename2plot <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == gene2plot, ]$external_gene_name
  colData(DDS)
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2) <- c("count", "condition")

  pdf(paste0(outdir, Project, "-DGE_", gene2plot,  "_", genename2plot,"_", TITLE, "_collated.pdf"),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({
    ggplot(t2, aes(x=condition, y=count, fill=condition)) +
      geom_violin(trim=TRUE, alpha=.5)  +
      geom_boxplot(width = 0.1, fill='white', outlier.shape=NA) +
      geom_point(position=position_jitter(w=0.1,h=0), alpha=0.5) +
      ggtitle(paste0(Project, " ::: ", gene2plot, " ::: ", genename2plot)) +
      xlab("") + ylab("Normalised count") +
      scale_fill_manual(name="Genotype", values = c("purple","purple4", "lightgreen", "green")) +
      theme(text = element_text(size=elementTextSize), legend.position="none") })
  dev.off()

  t2$samples   <- rownames(t2)
  colnames(t2) <- c("count", "condition", "samples")
  t2           <- t2[order(t2$condition),]
  t2$samples2 <- factor(t2$samples, as.character(t2$samples))

  pdf(paste0(outdir, "/", Project, "-DGE_", gene2plot,  "_", genename2plot, "_", TITLE, "_individual.pdf"),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ ggplot(t2, aes(x=samples2, y=count, fill=condition, group=condition)) +
      geom_bar(stat="identity", alpha=.5) +
      scale_fill_manual(name="Comparison", values = c("purple","purple4", "lightgreen", "green")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") +
      ggtitle(paste0(Project, " ::: ", gene2plot, " ::: ", genename2plot)) })
  dev.off()
  print(paste("Created plot for", gene2plot), sep=" ")

}
message("+------ Customised MA plot with/wo labeling some top DEGs and highlight key markers as suggested -----+")

functionCustomMAPlot <- function(results, Project, FigureID, Title, significance, l2fc, label_l2fc, xfavoritegenes, outdir) {

  rownames(results) <- make.names(results$ensembl_gene_id, unique=TRUE)
  labeldata.ann <- results[c(xfavoritegenes),]

  plt <- ggplot(data = results, aes(x=baseMean, y=log2FoldChange )) +
    geom_abline(intercept = l2fc, slope = 0, colour='red', alpha=0.25) +
    geom_abline(intercept = -l2fc, slope = 0, colour='red', alpha=0.25) +
    geom_point(size=0.5, alpha=0.5, col="black") +
    geom_point(data=subset(results, (padj <= significance & log2FoldChange >= l2fc)),  size=1, alpha=0.5,  col="red") +
    geom_point(data=subset(results, (padj <= significance & log2FoldChange <= -l2fc)), size=1, alpha=0.5,  col="blue") +
    geom_point(data=subset(labeldata.ann, log2FoldChange > 0), aes( x=baseMean, y=log2FoldChange), size=1.25, alpha=1.0, color='darkred',  shape=21, stroke=0.5) +
    geom_point(data=subset(labeldata.ann, log2FoldChange < 0), aes( x=baseMean, y=log2FoldChange), size=1.25, alpha=1.0, color='darkblue', shape=21, stroke=0.5) +
    geom_label_repel(data=subset(results, (padj <= significance & log2FoldChange >= label_l2fc)),
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name),
                     colour='black', point.padding = unit(0.25, "lines"),  size=3, segment.size = 1, segment.color = 'darkred',  nudge_x = 0, nudge_y=0) +
    geom_label_repel(data=subset(results, (padj <= significance & log2FoldChange <= -label_l2fc)),
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name),
                     colour='black', point.padding = unit(0.25, "lines"),  size=3, segment.size = 1, segment.color = 'darkblue',  nudge_x = 0, nudge_y=0) +
    geom_label_repel(data=labeldata.ann,
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name),
                     colour='purple', point.padding = unit(0.25, "lines"),  size=3, segment.size = 1, segment.color = 'purple',  nudge_x =1, nudge_y=-0.5) +
    scale_x_log10() +
    xlab("Mean Normalised Read Count") + ylab("log2 Fold Change") +
    ggtitle(paste(Project, " ", FigureID, "\n", Title, " [log2fc=", l2fc, ", sign=", significance, ", labels_log2fc=", label_l2fc, "]", sep=""))+
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


  pdf(paste0(outdir, "/", Project, "-MAplot_", FigureID, "_rm27_rm3_23_07_2019.pdf"),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ plt })
  dev.off()

  plt2 <- ggplot(data = results, aes(x=baseMean, y=log2FoldChange )) +
    geom_point(size=0.5, alpha=0.5, col="black") +
    geom_point(data=subset(results, (padj <= significance & log2FoldChange >= l2fc)),  size=1, alpha=0.5,  col="red") +
    geom_point(data=subset(results, (padj <= significance & log2FoldChange <= -l2fc)), size=1, alpha=0.5,  col="blue") +
    geom_point(data=subset(labeldata.ann, log2FoldChange > 0), aes( x=baseMean, y=log2FoldChange), size=1.25, alpha=1.0, color='darkred',  shape=21, stroke=0.5) +
    geom_point(data=subset(labeldata.ann, log2FoldChange < 0), aes( x=baseMean, y=log2FoldChange), size=1.25, alpha=1.0, color='darkblue', shape=21, stroke=0.5) +
    scale_x_log10() +
    xlab("Mean Normalised Read Count") + ylab("log2 Fold Change") +
    ggtitle(paste(Project, " ", FigureID, "\n", Title, " [log2fc=", l2fc, ", sign=", significance, "]", sep="")) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  pdf(paste0(outdir, "/", Project, "-MAplot_", FigureID, ".noannotations_rm27_rm3_23_07_2019.pdf"),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ plt2 })
  dev.off()

  return (plt)
}

message("+------                    Customised Volcano plot for different comparison groups -------------------+")
functionPlotDEVolcano <- function(results, sig_cut, logfc_cut, title,  xrange, yrange, topN, xlabel, ylabel, dnum, unum) {
  results       <- as.data.frame(results)
  results$genes <- results$external_gene_name
  results <- results[order(-results$log2FoldChange),]
  
  grob.down <- grobTree(textGrob(paste0("Down=", dnum), x=-3,  y=80, hjust=0,
                                 gp=gpar(col="blue", fontsize=13, fontface="bold")))
  grob.up<- grobTree(textGrob(paste0("Up=", dnum), x=3,  y=80, hjust=0,
                              gp=gpar(col="red", fontsize=13, fontface="bold")))
  # Plot
 
  volc.plt <- ggplot(data=results, aes(x=log2FoldChange, y=-log(padj), label=genes)) +
  
                    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
                    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
                    geom_hline(yintercept = -log(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
                    geom_point(data=results, alpha = 0.75, size=0.2, colour="grey") +
                    geom_point(data=subset(results, results$padj<=sig_cut & results$log2FoldChange >= logfc_cut ),
                               alpha = 0.75, size=0.8, colour="red") +
                    geom_point(data=subset(results, results$padj<=sig_cut & results$log2FoldChange <= -(logfc_cut)),
                               alpha = 0.75, size=0.8, colour="blue") +
   
                    geom_text_repel( data= subset(results, results$log2FoldChange > 1 & results$padj<=sig_cut & results$external_gene_name!="")[1:topN,],
                                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
                    geom_text_repel( data= tail(subset(results, results$log2FoldChange < (-1) & results$padj<=sig_cut & results$external_gene_name!=""), n=topN),
                                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
                    xlab(xlabel) + ylab(ylabel) +
                    scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) +
                    scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
                    theme(aspect.ratio=1) +
                    ggtitle(title) +
                    theme_update(plot.title = element_text(size=16, face="bold", hjust=0.5),
                               axis.title.x = element_text(size=12, face= "bold"),
                               axis.text.x = element_text(size=12, face="bold"),
                               axis.title.y.left = element_text(size=12, face= "bold"),
                               axis.text.y = element_text(size=12, face="bold")) +
                    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                  panel.background = element_rect(fill = "white", colour = NA)) +
    #annotation_custom(grob.down) +
    #annotation_custom(grob.up)
    
                    annotate(geom="text", x=-3, y=80, label=paste0("Down=", dnum), color="blue") +
                    annotate(geom="text", x=3, y=80, label=paste0("Up=", unum), color="red") 
                     
    return(volc.plt)
}

message("+-----------            Heatmap data sort and Heatmap plot function                     --------------+")

heatmap_data_fn <- function(res.file, psig.level=0.05, lfc.level=1, select.no, outfile){
  res.file1 <- res.file[!duplicated(res.file$ensembl_gene_id),];
  rownames(res.file1) <- res.file1$ensembl_gene_id
  res.flc <- rownames(head(res.file1[order(-abs(res.file1$log2FoldChange)),], n= select.no))
  res.sig.flc <- rownames(subset(res.file1, abs(res.file1$log2FoldChange) > lfc.level))

  res.padj <- rownames(head(res.file1[order(res.file1$padj),], n=select.no))
  res.sig.padj <- rownames(subset(res.file1, abs(res.file1$padj) < psig.level))
  res.flc.padj <- rownames(subset(res.file1, abs(res.file1$log2FoldChange) > lfc.level & res.file1$padj < psig.level))
  save(res.flc, res.sig.flc, res.padj, res.sig.padj, res.flc.padj, file = outfile)
}
##
heatmap_mat_fn <- function(rlds, genes2plot, ensE, samTable, del.cols, sel.cols, outfile, width, height){
  mat.rows <- match(genes2plot, row.names(rlds))
  mat.rows <- mat.rows[!is.na(mat.rows)]
  mat1 <- assay(rlds)[mat.rows,]
  mat1.df <- data.frame(ensembl_gene_id=rownames(mat1), mat1)
  mat1.ann <- merge(mat1.df, ensE, by ="ensembl_gene_id")
  mat1.ann$dups <- duplicated(mat1.ann$ensembl_gene_id)
  mat1.ann.dedup <- subset(mat1.ann, mat1.ann$dups == F)
  rownames(mat1.ann.dedup) <- paste0(mat1.ann.dedup$ensembl_gene_id, "::", mat1.ann.dedup$external_gene_name)
  # rownames(mat1.ann.dedup) <-  mat1.ann.dedup$external_gene_name
  mat2 <- mat1.ann.dedup[,-c(del.cols)]
  samTable <- as.data.frame(samTable)
  samTable <- samTable[order(samTable$Labels),]
  rownames(samTable) <- samTable$sampleName
  colnames(mat2) <- samTable$sampleName

  annotation_col <- as.data.frame(samTable[, sel.cols])
  rownames(annotation_col) <- samTable$sampleName
  Treatment <- annotation_col$condition
  ann_colors <- list(Treatment = c(TX="steelblue3", SHAM="darkred" ))
  pdf(outfile, onefile=FALSE, width=width, height=height)
  par(bg=NA)
  pheatmap(mat2,  annotation_col = annotation_col, annotation_colors = ann_colors,
           fontsize=7, fontsize_row=2, show_rownames=F, cluster_cols = T,
           cluster_rows=T, cutree_cols = 2)
  dev.off()

}
## merge treatment and age together to plot a heatmap.
heatmap_mat_fn_age <- function(rlds, genes2plot, ensE, samTable, outfile, width, height){
  mat.rows <- match(genes2plot, row.names(rlds))
  mat.rows <- mat.rows[!is.na(mat.rows)]
  mat1 <- assay(rlds)[mat.rows,]
  mat1.df <- data.frame(ensembl_gene_id=rownames(mat1), mat1)
  mat1.ann <- merge(mat1.df, ensE, by ="ensembl_gene_id")
  mat1.ann$dups <- duplicated(mat1.ann$ensembl_gene_id)
  mat1.ann.dedup <- subset(mat1.ann, mat1.ann$dups == F)
  rownames(mat1.ann.dedup) <- paste0(mat1.ann.dedup$ensembl_gene_id, "::", mat1.ann.dedup$external_gene_name)
  keep.cols <- grep("^s_*", colnames(mat1.ann.dedup), value = F)
  mat2 <- mat1.ann.dedup[,keep.cols]

  samTable <- samTable[order(samTable$Labels),]
  colnames(mat2) <- samTable$sample_id
  sel.cols <- which(colnames(samTable)=="condition" | colnames(samTable)=="GestationalAge")


  samTable <- as.data.frame(samTable)
  samTable <- samTable[order(samTable$Labels),]
  rownames(samTable) <- samTable$sample_id
  colnames(mat2) <- samTable$sample_id

  annotation_col <- as.data.frame(samTable[,which(colnames(samTable)=="CombineCondition")])
  rownames(annotation_col) <- samTable$sample_id
  colnames(annotation_col) <- "CombineCondition"
  Treatment <- annotation_col$CombineCondition
  ann_colors <- list(Treatment = c(TX_129="steelblue3", TX_143="darkgreen", SHAM_129="orange",SHAM_143="darkred" ))
  pdf(outfile, onefile=FALSE, width=width, height=height)
  par(bg=NA)
  pheatmap(mat2,  annotation_col = annotation_col, annotation_colors = ann_colors,
           fontsize=7, fontsize_row=7, show_rownames=F, cluster_cols = T, cutree_cols = 2,
           cluster_rows=T, treeheight_row = 0)
  dev.off()

}

message("+-------------------------------------------------------------------------------+")
message("+--  Differential of Differentials, compare 143 TX vs SHAM to 129 TX vs SHAM  --+")
message("+-------------------------------------------------------------------------------+")

functionPlotDECorrelation <- function(clustA, clustB, clustA.lab, clustB.lab, sig_cut, logfc_cut, topN, out.dir, nameTB, nameTG, nameTP, nameTY) {

  rownames(clustA)  <- make.names(clustA$ensembl_gene_id, unique=T)
  rownames(clustB)  <- make.names(clustB$ensembl_gene_id, unique=T)

  clustA <- clustA[,c("log2FoldChange",  "padj",  "external_gene_name")]
  clustB <- clustB[,c("log2FoldChange",  "padj",  "external_gene_name")]

  colnames(clustA)  <- c("log2FoldChange.A",  "padj.A",  "external_gene_name.A")
  colnames(clustB)  <- c("log2FoldChange.B",  "padj.B",  "external_gene_name.B")

  clustA.topN <- subset(clustA, clustA$log2FoldChange.A > 0      & clustA$padj.A<=sig_cut)[1:topN,]
  clustA.botN <- tail(subset(clustA, clustA$log2FoldChange.A < 0 & clustA$padj.A<=sig_cut),topN)

  clustB.topN <- subset(clustB, clustB$log2FoldChange.B > 0      & clustB$padj.B<=sig_cut)[1:topN,]
  clustB.botN <- tail(subset(clustB, clustB$log2FoldChange.B < 0 & clustB$padj.B<=sig_cut),topN)

  compare.A.B <- merge(clustA, clustB, by='row.names', all=TRUE)
  compare.A.B[ is.na(compare.A.B)] <- 0

  compare.A.B$colour[(compare.A.B$padj.A <= sig_cut & compare.A.B$padj.B <= sig_cut)] <- "purple"
  # & abs(compare.A.B$log2FoldChange.A) >= logfc_cut & abs(compare.A.B$log2FoldChange.B) >= logfc_cut)] <- "purple"

  compare.A.B$colour[(compare.A.B$padj.A > sig_cut  & compare.A.B$padj.B <= sig_cut)] <- "blue"
  # & abs(compare.A.B$log2FoldChange.A) < logfc_cut & abs(compare.A.B$log2FoldChange.B) >= logfc_cut)] <- "blue"
  compare.A.B$colour[(compare.A.B$padj.A == 0       & compare.A.B$padj.B <= sig_cut)] <- "blue"
  # & abs(compare.A.B$log2FoldChange.A) < logfc_cut & abs(compare.A.B$log2FoldChange.B) >= logfc_cut)] <- "blue"

  compare.A.B$colour[(compare.A.B$padj.A <= sig_cut & compare.A.B$padj.B > sig_cut)] <- "darkgreen"
  # & abs(compare.A.B$log2FoldChange.A) >= logfc_cut & abs(compare.A.B$log2FoldChange.B) < logfc_cut)] <- "darkgreen"
  compare.A.B$colour[(compare.A.B$padj.A <= sig_cut & compare.A.B$padj.B == 0      )] <- "darkgreen"
  #& abs(compare.A.B$log2FoldChange.A) >= logfc_cut & abs(compare.A.B$log2FoldChange.B) < logfc_cut)] <- "darkgreen"

  compare.A.B$colour[(compare.A.B$padj.A > sig_cut & compare.A.B$padj.B > sig_cut)] <- "grey"
  compare.A.B$colour[(compare.A.B$padj.A > sig_cut & compare.A.B$padj.B == 0 )]     <- "grey"
  compare.A.B$colour[(compare.A.B$padj.A == 0      & compare.A.B$padj.B > sig_cut)] <- "grey"
  compare.A.B$colour[(compare.A.B$padj.A == 0      & compare.A.B$padj.B == 0 )] <- "grey"


  compare.A.B$colour[( abs(compare.A.B$log2FoldChange.A) < logfc_cut & abs(compare.A.B$log2FoldChange.B) < logfc_cut)] <- "grey"

  clustB.blue.all      <- subset(compare.A.B, compare.A.B$colour == "blue" & compare.A.B$padj.B < sig_cut & compare.A.B$padj.B!=0)
  clustB.blue.all.num <- dim(clustB.blue.all)[1]
  write.csv(clustB.blue.all, file=paste0(out.dir, "/", Project, "-", nameTB, clustB.blue.all.num, "_sigP", sig_cut, "_l2FC", logfc_cut, "_all.ann.csv"))


  clustB.blue      <- subset(compare.A.B, colour == "blue" & padj.B < sig_cut & padj.B!=0 & abs(log2FoldChange.B) >= logfc_cut)
  clustB.blue.topN <- clustB.blue[order(clustB.blue$log2FoldChange.B,decreasing=TRUE),]
  clustB.blue.topN <- clustB.blue.topN[c(1:5),]
  clustB.blue.botN <- clustB.blue[order(clustB.blue$log2FoldChange.B,decreasing=FALSE),]
  clustB.blue.botN <- clustB.blue.botN[c(1:5),]
  clustB.blue.num <- dim(clustB.blue)[1]
  write.csv(clustB.blue, file=paste0(out.dir, "/", Project, "-", nameTB, clustB.blue.num, "_sigP", sig_cut, "_l2FC", logfc_cut, "_sig.ann.csv"))


  clustA.green.all      <- subset(compare.A.B, compare.A.B$colour == "darkgreen" & compare.A.B$padj.A < sig_cut& compare.A.B$padj.A!=0)
  clustA.green.all.num <- dim(clustA.green.all)[1]
  write.csv(clustA.green.all, file=paste0(out.dir, "/", Project, "-", nameTG,clustA.green.all.num, "_sigP", sig_cut, "_l2FC", logfc_cut, "_all.ann.csv"))


  clustA.green      <- subset(compare.A.B, colour == "darkgreen" & padj.A < sig_cut & padj.A!=0 & abs(log2FoldChange.A) >= logfc_cut)
  clustA.green.topN <- clustA.green[order(clustA.green$log2FoldChange.A,decreasing=TRUE),]
  clustA.green.topN <- clustA.green.topN[c(1:5),]
  clustA.green.botN <- clustA.green[order(clustA.green$log2FoldChange.A,decreasing=FALSE),]
  clustA.green.botN <- clustA.green.botN[c(1:5),]
  clustA.green.num <- dim(clustA.green)[1]
  write.csv(clustA.green, file=paste0(out.dir, "/", Project, "-",nameTG, clustA.green.num, "_sigP", sig_cut, "_l2FC", logfc_cut, "_sig.ann.csv"))

  clustA.B.purple.all <- subset(compare.A.B, colour == "purple")
  clustA.B.purple.all.num <- dim(clustA.B.purple.all)[1]
  write.csv(clustA.B.purple.all, file=paste0(out.dir, "/", Project, "-", nameTP, clustA.B.purple.all.num, "_sigP", sig_cut, "_l2FC", logfc_cut, "_all.ann.csv"))

  clustA.B.purple <- subset(compare.A.B, colour == "purple" & padj.A < sig_cut & padj.B < sig_cut & abs(log2FoldChange.A) >= logfc_cut & abs(log2FoldChange.B)>=logfc_cut)
  clustA.B.purple.num <- dim(clustA.B.purple)[1]
  write.csv(clustA.B.purple, file=paste0(out.dir, "/", Project, "-", nameTP, clustA.B.purple.num, "_sigP", sig_cut, "_l2FC", logfc_cut, "_sig.ann.csv"))
  clustA.B.grey <- subset(compare.A.B, colour == "grey")
  clustA.B.grey.num <- dim(clustA.B.grey)[1]
  write.csv(clustA.B.grey, file=paste0(out.dir, "/", Project, "-", nameTY, clustA.B.grey.num, "_sigP", sig_cut, "_l2FC", logfc_cut, "_all.ann.csv"))


  compare.A.B$label <- 0
  compare.A.B[compare.A.B$Row.names %in% unlist(rownames(clustA.topN)), ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(rownames(clustA.botN)), ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(rownames(clustB.topN)), ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(rownames(clustB.botN)), ]$label  <- 1

  compare.A.B[compare.A.B$Row.names %in% unlist(clustB.blue.topN$Row.names),  ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(clustB.blue.botN$Row.names),  ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(clustA.green.topN$Row.names), ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(clustA.green.botN$Row.names), ]$label  <- 1

  #minFC <- round(min(compare.A.B$log2FoldChange.A, compare.A.B$log2FoldChange.B))
  #maxFC <- ceiling(max(compare.A.B$log2FoldChange.A, compare.A.B$log2FoldChange.B))
  minFC <- -8 #min( compare.A.B$log2FoldChange.A, compare.A.B$log2FoldChange.B  )
  maxFC <-  8 #max( compare.A.B$log2FoldChange.A, compare.A.B$log2FoldChange.B  )

  cor.plt <- ggplot(data=compare.A.B, aes(x=log2FoldChange.A, y=log2FoldChange.B, colour=colour,
                                          label=external_gene_name.A)) +
    geom_vline(xintercept = -(logfc_cut), linetype="dashed", colour="black", size=0.25, alpha=0.5) +
    geom_vline(xintercept = logfc_cut,    linetype="dashed", colour="black", size=0.25, alpha=0.5) +
    geom_hline(yintercept = -(logfc_cut), linetype="dashed", colour="black", size=0.25, alpha=0.5) +
    geom_hline(yintercept = logfc_cut,    linetype="dashed", colour="black", size=0.25, alpha=0.5) +
    geom_abline(intercept = 0, linetype="solid", colour="black", alpha=0.5) +
    geom_point(alpha=0.5, size=0.5) +
    geom_label_repel( data=subset(compare.A.B, compare.A.B$label == 1 ), force=10, fill='white',
                      show.legend = FALSE, nudge_x=0.0, nudge_y=0.0, segment.size = 0.25, size=3) +
    coord_fixed() +
    scale_x_continuous(limits=c(minFC, maxFC), breaks=seq(minFC, maxFC,2)) +
    scale_y_continuous(limits=c(minFC, maxFC), breaks=seq(minFC, maxFC,2)) +
    xlab(paste0("log2FC (", clustA.lab, ") ")) +
    ylab(paste0("log2FC (", clustB.lab, ") ")) +
    scale_colour_manual(name="", values=c("purple"="purple", "blue"="blue",
                                          "darkgreen"="darkgreen", "grey"="grey"),
                        labels=c("purple"=paste0(clustA.lab, " & ", clustB.lab),
                                 "blue"=paste0(clustB.lab, " only"),
                                 "darkgreen"=paste0(clustA.lab, " only"),
                                 "grey" = "nonSig")) +
    ggtitle(paste0(clustA.lab," Vs ", clustB.lab, " [log2FC=", logfc_cut,"]")) +
    theme(text=element_text(family="sans"), legend.position="bottom", aspect.ratio=1) +
    theme_classic()


  return(cor.plt)
}
