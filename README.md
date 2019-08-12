# Thyroid deficiency in utero stimulates white adipose tissue growth and impairs brown adipose thermogenic capacity

Protemics analysis for Mouse Placenta Data from Scaffold results. (P147_2_100718.sf3) Five samples have been analysed, 1467 proteins have been identified based on the scaffold default settings, protein threshold 95%, min number of peptides is 2 and peptide threshold is 95%. To have a high
quality data set for annotation homolog analysis, we choose the proteins which have been
identified in at least 4 samples. We only consider the mouse protein identification and find the homolog protein in human. There are 921 proteins in the list has been selected initially based on the conditions above. </br>

Our purpose is to identify the homolog proteins between mouse and human, we would like to download the data sets which relating to mouse/human placenta. And perform the expression analysis, then make a general normalisation method to make these data sets comparable. If one gene expressed evenly in both human and mouse data sets, then this protein for Normal Placentas is a good candidate, otherwise, the gene may correspond to some diseases. </br>

Later in May, we had FACS peptide protein data, which has pellet (particles can sediment at the bottom of the tube and this isolated specimen) and supernatant (the remaining solution which can be further processed or analysed.) There are three samples and the data format with peptides counts. (S3/P3, S4/P4, S5/P5) In order to have a good understanding the previous scaffold data protocol and the FACs method, we would like to merge supernatant and pellet together and do the same filtering as scaffold data above. </br>


## Publication

## Abstract

## Data Processing
Data were aligned to sheep genome (Ovis_aries.Oar_v3.1) with STAR (2.5.1b_modified) with --single option. Alignments and QC were processed using custom ClusterFlow (v0.5dev) pipelines and assessed using MultiQC(0.9dev).
Gene quantification was determined with Subread (1.5.0-p2) function featureCounts function. Differential gene expression was performed with DESeq2 package (1.22.2, R 3.5.3) and with the same package read counts were normalised
on the estimated size factors. </bt>

Gene ontology (GO) and KEGG pathway analysis were performed using clusterProfiler (3.10.1).  </bt>

<table>
  <tr>
    <th>Resource</th><th>URL</th>
  </tr>
  <tr>
    <td>Oar_v3.1</td><td>[Link](https://www.ensembl.org/Ovis_aries/Info/Index)</td>
  </tr>
  <tr>
    <td>FastQC</td><td>[Link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)</td>
  </tr>
  <tr>
    <td>Trim_galore</td><td>[Link](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)</td>
  </tr>
  <tr>
    <td>STAR</td><td>[DOI](https://academic.oup.com/bioinformatics/article/29/1/15/272537)</td>
  </tr>
  <tr>
    <td>Feature_counts</td><td>[DOI](https://academic.oup.com/bioinformatics/article/30/7/923/232889)</td>
  </tr>
  <tr>
    <td>ClusterFlow</td><td>[DOI](https://f1000research.com/articles/5-2824/v2)</td>
  </tr>
  <tr>
    <td>MultiQC</td><td>[DOI](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)</td>
  </tr>
</table>



### Stage 1: Normal Placentas

  - [x]L1: Comparing to mouse Placenta (RNA-Seq, microarray)
  - [x]L2: Comparing to human Placenta (RNA-Seq, microarray)
  - [?]L3: Comparing to Lz and Jz parts of the rat placenta RNA-Seq, the data is not available for the analysis at this moment, only one reference from stanford but not useful;
  - [?]L4: (\*) Comparing to specific cell type in the human placenta (Syncytiotrophoblast,
    Amanda mentioned that it is easy in human but may not easy for mouse, need to check later)
  - [?]L5: Comparing to RNA-Seq from secreted placental exosome
  - [?]Single-cell RNAseq data set check </br>
ps: The analysis we perform the comparison with public data for normal
### Stage 2 : Disease Placentas (Pregnancy conditioners)

  - [x] L1: Disease public data sets gathering; conditions are PEs(Preeclampsia), GDM(gestational diabetes), SGA(small gestational age), LGA(large gestational age), and recurrent miscarriage;
  - [x] L2: Overlap genes expression across all available data sets to check the expressed genes for all of the possible complications;

### Stage 3: Igf2 Control vs Igf2 Knock out (CTR_ans48_0005--Igf2)

- [x] Default settings from scaffold for the proteins selection are, protein threshold: 95%, minimum number of peptides are 2 and peptide threshold is 95%. (1419 proteins are selected after this filtering to do our next analysis)
  - [x] Select the subset of Protein data sets, which includes
    - [x] existing in both condition, appeared at least in two replicates for each condition;
    - [x] Or existing in either condition, have at least two replicates, these proteins may specifically relating to the conditions.
  - [x] Fisher-g-test for all of the select subsets proteins on the foldchange, to detect the up/down regulated proteins and then linked with genes;
    - [?] The replicates spectra counts are not uniform within the conditions, do we need to consider this abnormal distribution of counts or just ignore it using the mean value?
    - [?] If using the mean values, for some proteins, they may have unequal number of replicates in both conditions, e.g. Control 2 samples, and KO 3 samples;
    - [?] How to put the test also having the consideration of peptides counts?

## Scripts
## Tables&Figures
