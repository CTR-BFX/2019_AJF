# Thyroid deficiency in utero stimulates white adipose tissue growth and impairs brown adipose thermogenic capacity




## Publication

## Abstract

## Data Processing
Data were aligned to sheep genome (Oar_v3.1) with STAR (2.5.1b_modified) with --single option. Alignments and QC were processed using custom ClusterFlow (v0.5dev) pipelines and assessed using MultiQC(0.9dev).
Gene quantification was determined with Subread (1.5.0-p2) function featureCounts function. Differential gene expression was performed with DESeq2 package (1.22.2, R 3.5.3) and with the same package read counts were normalised
on the estimated size factors. </bt>

Gene ontology (GO) and KEGG pathway analysis were performed using clusterProfiler (3.10.1).  </bt>

Resource       | URL
-------------- | --------------
Oar_v3.1       | [Link](https://www.ensembl.org/Ovis_aries/Info/Index)
FastQC         | [Link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
Trim_galore    | [Link](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
STAR           | [DOI](https://academic.oup.com/bioinformatics/article/29/1/15/272537)
Feature_Counts | [DOI](http://dx.doi.org/10.1093/bioinformatics/btt656)
ClusterFlow    | [DOI](http://dx.doi.org/10.12688/f1000research.10335.2)
MultiQC        | [DOI](http://dx.doi.org/10.1093/bioinformatics/btw354)


## Scripts to reproduce paper figures and Tables

All files are provides in this repository with the exception of the GTF file for the sheep reference genome. The GTF file can be downloaded from https://www.ensembl.org/Ovis_aries/Info/Index. </bt>

The bash script for align and QC control is named as below.

                                       ClusterFlow_FTSEQ.sh

All the featureCounts files are in the Feature_Counts folder. These scripts can be run interactively in R-studio or as a batch using **Rscripts** by changing the working directory. The *Backstge_Functions_19_07_2019.R* is the useful functions which will be called in the main script *DESeq2_Analysis.R*. </bt>

### Additional Methods Information

## Session Information

Details for the R version and the corresponding used packages to creat all Figures and Tables is given below.

````
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.14.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] ComplexHeatmap_2.1.0        UpSetR_1.4.0               
 [3] ggalt_0.4.0                 ggthemes_4.2.0             
 [5] calibrate_1.7.2             MASS_7.3-51.4              
 [7] dplyr_0.8.3                 plyr_1.8.4                 
 [9] biomaRt_2.38.0              reshape2_1.4.3             
[11] ggrepel_0.8.1               pheatmap_1.0.12            
[13] cowplot_1.0.0               gplots_3.0.1.1             
[15] RColorBrewer_1.1-2          ggplot2_3.2.0              
[17] DESeq2_1.22.2               SummarizedExperiment_1.12.0
[19] DelayedArray_0.8.0          BiocParallel_1.16.6        
[21] matrixStats_0.54.0          ensembldb_2.6.8            
[23] AnnotationFilter_1.6.0      GenomicFeatures_1.34.8     
[25] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
[27] circlize_0.4.6              AnnotationHub_2.14.5       
[29] pathview_1.22.3             org.Hs.eg.db_3.7.0         
[31] enrichplot_1.2.0            org.Bt.eg.db_3.7.0         
[33] rgdal_1.4-4                 sp_1.3-1                   
[35] UniProt.ws_2.22.0           RCurl_1.95-4.12            
[37] bitops_1.0-6                RSQLite_2.1.2              
[39] GEOquery_2.50.5             GSEABase_1.44.0            
[41] graph_1.60.0                annotate_1.60.1            
[43] XML_3.98-1.20               AnnotationDbi_1.44.0       
[45] IRanges_2.16.0              S4Vectors_0.20.1           
[47] Biobase_2.42.0              BiocGenerics_0.28.0        
[49] DOSE_3.8.2                  clusterProfiler_3.10.1     

loaded via a namespace (and not attached):
  [1] tidyselect_0.2.5              htmlwidgets_1.3              
  [3] munsell_0.5.0                 withr_2.1.2                  
  [5] colorspace_1.4-1              GOSemSim_2.8.0               
  [7] knitr_1.23                    rstudioapi_0.10              
  [9] Rttf2pt1_1.3.7                KEGGgraph_1.42.0             
 [11] urltools_1.7.3                GenomeInfoDbData_1.2.0       
 [13] polyclip_1.10-0               bit64_0.9-7                  
 [15] farver_1.1.0                  vctrs_0.2.0                  
 [17] xfun_0.8                      BiocFileCache_1.6.0          
 [19] R6_2.4.0                      clue_0.3-57                  
 [21] locfit_1.5-9.1                fgsea_1.8.0                  
 [23] gridGraphics_0.4-1            assertthat_0.2.1             
 [25] promises_1.0.1                scales_1.0.0                 
 [27] ggraph_1.0.2                  nnet_7.3-12                  
 [29] gtable_0.3.0                  ash_1.0-15                   
 [31] rlang_0.4.0                   zeallot_0.1.0                
 [33] genefilter_1.64.0             GlobalOptions_0.1.0          
 [35] splines_3.5.1                 extrafontdb_1.0              
 [37] rtracklayer_1.42.2            lazyeval_0.2.2               
 [39] acepack_1.4.1                 europepmc_0.3                
 [41] checkmate_1.9.4               BiocManager_1.30.4           
 [43] yaml_2.2.0                    backports_1.1.4              
 [45] httpuv_1.5.1                  qvalue_2.14.1                
 [47] Hmisc_4.2-0                   extrafont_0.17               
 [49] tools_3.5.1                   ggplotify_0.0.3              
 [51] ggridges_0.5.1                Rcpp_1.0.2                   
 [53] base64enc_0.1-3               progress_1.2.2               
 [55] zlibbioc_1.28.0               purrr_0.3.2                  
 [57] prettyunits_1.0.2             rpart_4.1-15                 
 [59] GetoptLong_0.1.7              viridis_0.5.1                
 [61] cluster_2.1.0                 magrittr_1.5                 
 [63] data.table_1.12.2             DO.db_2.9                    
 [65] triebeard_0.3.0               ProtGenerics_1.14.0          
 [67] hms_0.5.0                     mime_0.7                     
 [69] xtable_1.8-4                  gridExtra_2.3                
 [71] shape_1.4.4                   compiler_3.5.1               
 [73] maps_3.3.0                    tibble_2.1.3                 
 [75] KernSmooth_2.23-15            crayon_1.3.4                 
 [77] htmltools_0.3.6               later_0.8.0                  
 [79] Formula_1.2-3                 tidyr_0.8.3                  
 [81] geneplotter_1.60.0            DBI_1.0.0                    
 [83] tweenr_1.0.1                  proj4_1.0-8                  
 [85] dbplyr_1.4.2                  rappdirs_0.3.1               
 [87] Matrix_1.2-17                 readr_1.3.1                  
 [89] gdata_2.18.0                  igraph_1.2.4.1               
 [91] pkgconfig_2.0.2               rvcheck_0.1.3                
 [93] GenomicAlignments_1.18.1      foreign_0.8-72               
 [95] xml2_1.2.1                    XVector_0.22.0               
 [97] stringr_1.4.0                 digest_0.6.20                
 [99] Biostrings_2.50.2             fastmatch_1.1-0              
[101] htmlTable_1.13.1              curl_4.0                     
[103] shiny_1.3.2                   Rsamtools_1.34.1             
[105] gtools_3.8.1                  rjson_0.2.20                 
[107] jsonlite_1.6                  viridisLite_0.3.0            
[109] limma_3.38.3                  pillar_1.4.2                 
[111] lattice_0.20-38               KEGGREST_1.22.0              
[113] httr_1.4.0                    survival_2.44-1.1            
[115] GO.db_3.7.0                   interactiveDisplayBase_1.20.0
[117] glue_1.3.1                    png_0.1-7                    
[119] bit_1.1-14                    Rgraphviz_2.26.0             
[121] ggforce_0.2.2                 stringi_1.4.3                
[123] blob_1.2.0                    latticeExtra_0.6-28          
[125] caTools_1.17.1.2              memoise_1.1.0    

````

## SampleTable

                                      SampleTable.csv
## Links
Description   | URL
------------- | ----------
Publication   | [[Developmental Biology, <b>x</b>, xxxx]](https://www.xxxx.com/articles/xxxx) [[DOI]](https://doi.org/xxxx) <br> [[bioRxiv]](https://www.biorxiv.org/content/early/2018/05/24/330068) [[DOI]](https://doi.org/10.1101/330068)
Raw Data      | ArrayExpress EMBL-EBI [E-MTAB-xxxx](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-xxxx)
CTR Bioinformatics | [CTR-BFX](https://www.trophoblast.cam.ac.uk/Resources/BioInformatics)

## Contact
Contact Xiaohui Zhao (xz289@cam.ac.uk) for bioinformatics related queries.
