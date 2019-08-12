#!/usr/local/bin/bash
#
# CTR_ajf1005_0001 ::: RNAseq clusterflow+htseq_count
#
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
# initial install multiqc, and HTSeq
# 
# pip install multiqc
# pip install HTSeq

#--------------------------------------------------------#
#                   clusterflow                          #
#--------------------------------------------------------#
ajfdir=/storage/CTR-Projects/CTR_ajf1005/CTR_ajf1005_0001
cd $ajfdir/OriData_QC
module load clusterflow/v0.5dev
cf --genome Oar_v3.1 --single fastq_star *.fq.gz

#--------------------------------------------------------#
#                   multiqc report                       #
#           cluster folder for packdir                   #                 
#--------------------------------------------------------#
packdir=/usr/local/bin
cd $ajfdir
$packdir/multiqc -f -c CTR_ajf1005_0001.config.yaml  OriData_QC/

#--------------------------------------------------------#
#   htseq-count, clusterflow version                     #
#--------------------------------------------------------#
inFiles=$ajfdir/OriData_QC/*.bam
cf --genome Oar_v3.1 htseq_counts $ajfdir/OriData_QC/*.bam
 

  
