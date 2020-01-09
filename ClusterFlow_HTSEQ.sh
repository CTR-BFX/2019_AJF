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
 

#--------------------------------------------------------#
# mergeBam files and sort, index by condition and age    #
#--------------------------------------------------------#

cd OriData_QC/
samtools merge ./MergeBam/SHAM_143_merge.bam \
                SLX-14347.D701_D504.*.star.bam \
                SLX-14347.D701_D507.*.star.bam \
                SLX-14347.D702_D504.*.star.bam \
                SLX-14347.D702_D505.*.star.bam \
                SLX-14347.D705_D506.*.star.bam \
                SLX-14347.D707_D505.*.star.bam \
                SLX-14347.D708_D506.*.star.bam \
                SLX-14347.D712_D504.*.star.bam \
                SLX-14347.D712_D506.*.star.bam 
             
samtools merge ./MergeBam/SHAM_129_merge.bam \
                SLX-14347.D701_D505.*.star.bam \
                SLX-14347.D703_D504.*.star.bam \
                SLX-14347.D703_D506.*.star.bam \
                SLX-14347.D704_D506.*.star.bam \
                SLX-14347.D706_D506.*.star.bam \
                SLX-14347.D709_D504.*.star.bam \
                SLX-14347.D709_D506.*.star.bam \
                SLX-14347.D710_D504.*.star.bam \
                SLX-14347.D711_D504.*.star.bam 

samtools merge ./MergeBam/TX_143_merge.bam \
                SLX-14347.D701_D506.*.star.bam \
                SLX-14347.D704_D505.*.star.bam \
                SLX-14347.D705_D505.*.star.bam \
                SLX-14347.D706_D504.*.star.bam \
                SLX-14347.D706_D505.*.star.bam \
                SLX-14347.D707_D506.*.star.bam \
                SLX-14347.D708_D505.*.star.bam \
                SLX-14347.D709_D505.*.star.bam \
                SLX-14347.D710_D506.*.star.bam \
                SLX-14347.D711_D505.*.star.bam 

samtools merge ./MergeBam/TX_129_merge.bam \
                SLX-14347.D702_D506.*.star.bam \
                SLX-14347.D703_D505.*.star.bam \
                SLX-14347.D704_D504.*.star.bam \
                SLX-14347.D705_D504.*.star.bam \
                SLX-14347.D707_D504.*.star.bam \
                SLX-14347.D708_D504.*.star.bam \
                SLX-14347.D710_D505.*.star.bam \
                SLX-14347.D711_D506.*.star.bam \
                SLX-14347.D712_D505.*.star.bam

for names in SHAM_143 SHAM_129 TX_143 TX_129;
do
  samtools sort ${names}_merge.bam -o ${names}_merge.srt.bam
  samtools index ${names}_merge.srt.bam
done
 
#--------------------------------------------------------#
#   featureCounts for THRA 1, THRA 2                     #
#--------------------------------------------------------#

cd ./MergeBam

for file in SHAM_129 SHAM_143 TX_129 TX_143
  do
  echo ${file}_merge.srt.bam
  featureCounts -T 1 -t exon -g gene_id -a THRA_1.gtf -o ${file}.THRA_1.fcounts ${file}_merge.srt.bam
  featureCounts -T 1 -t exon -g gene_id -a THRA_2.gtf -o ${file}.THRA_2.fcounts ${file}_merge.srt.bam
done


