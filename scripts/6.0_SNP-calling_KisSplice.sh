#!/bin/bash


#Rna-seq using STAR
#http://kissplice.prabi.fr/TWAS/
#./6.0_SNP-calling.sh > 6.0_SNP-calling.sh.log 2>&1


echo "run script for rna-seq-analysis"
##########################################
##########################################
INDIR=/home/gala0002/proj/proj_dir/NG-14833_1.1_sort-trim/
work_dir1=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Tit
work_dir2=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Reg
work_dir3=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Pas

out_dir1=/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_Tit/
out_dir2=/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_Reg/
out_dir3=/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_Pas/


out_dir1_TD=/home/gala0002/proj/proj_dir/NG-14833_2.0_TD_Tit
out_dir2_TD=/home/gala0002/proj/proj_dir/NG-14833_2.0_TD_Reg
out_dir3_TD=/home/gala0002/proj/proj_dir/NG-14833_2.0_TD_Pas

Trinity="/bioinfo/trinityrnaseq-Trinity-v2.8.5/Trinity"
kissplice="/bioinfo/KISSPLICE/kissplice-2.4.0-p1/bin/kissplice"
kissplice2reftrans="/bioinfo/KISSPLICE/kissplice2reftranscriptome-1.3.3/kissplice2reftranscriptome"
TransDecoder="/bioinfo/TransDecoder-TransDecoder-v5.5.0/"
TMP="/home/gala0002/proj/proj_dir/tmp/"

###########################################
#De novo assembly of reads using Trinity

#NG-14833_Tit

# nice -n 5 ${Trinity} --seqType fq --SS_lib_type RF \
# --left ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_1.fq.gz \
# --right ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_2.fq.gz \
# --CPU 25 --max_memory 20G --output ${work_dir1} --workdir ${work_dir1}

#NG-14833_Reg

# nice -n 5 ${Trinity} --seqType fq --SS_lib_type RF \
# --left ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_1.fq.gz \
# --right ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_2.fq.gz \
# --CPU 25 --max_memory 20G --output ${work_dir2} --workdir ${work_dir2}

#NG-14833_Pas

# nice -n 5 ${Trinity} --seqType fq --SS_lib_type RF \
# --left ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_1.fq.gz \
# --right ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_2.fq.gz \
# --CPU 30 --max_memory 50G --output ${work_dir3} --workdir ${work_dir3}

###########################################
#Running KisSplice

#${Kissplice} -s 1 -k 41 --experimental \
#-r cond1_replicat1_R1.fastq -r cond1_replicat1_R2.fastq -r cond1_replicat2_R1.fastq -r cond1_replicat2_R2.fastq \
#-r cond2_replicat1_R1.fastq -r cond2_replicat1_R2.fastq -r cond2_replicat2_R1.fastq -r cond2_replicat2_R2.fastq \


#NG-14833_Tit

# ulimit -s unlimited
# nice -n 5 ${kissplice} -s 1 -k 41 --experimental \
# -r ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2.fq.gz -r ${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_1.fq.gz -r ${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_2.fq.gz -r ${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_2.fq.gz \
# -r ${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_2.fq.gz -r ${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_2.fq.gz -r ${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_1.fq.gz -r ${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_2.fq.gz \
# -o ${out_dir1} -d ${TMP} -t 25


#NG-14833_Reg
# cd $out_dir2

# ulimit -s unlimited
# nice -n 5 ${kissplice} -s 1 -k 41 --experimental \
# -r ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2.fq.gz -r ${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_1.fq.gz -r ${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_2.fq.gz -r ${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_1.fq.gz -r ${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_2.fq.gz \
# -r ${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_1.fq.gz -r ${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_2.fq.gz -r ${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_2.fq.gz -r ${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_2.fq.gz \
# -o ${out_dir2} -d ${TMP} -t 20


# #NG-14833_Pas
# cd $out_dir3

# ulimit -s unlimited
# nice -n 5 ${kissplice} -s 1 -k 41 --experimental \
# -r ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2.fq.gz -r ${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_2.fq.gz -r ${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_2.fq.gz \
# -r ${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_2.fq.gz -r ${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_2.fq.gz -r ${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_1.fq.gz -r ${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_2.fq.gz \
# -o ${out_dir3} -d ${TMP} -t 20


###########################################
#Running Transdecoder

#Transdecoder.LongOrfs -t Trinity.fasta
#TransDecoder.Predict -t Trinity.fasta

# #NG-14833_Tit
# ${TransDecoder}TransDecoder.LongOrfs -t ${work_dir1}/Trinity.fasta 
# ${TransDecoder}TransDecoder.Predict -t ${work_dir1}/Trinity.fasta 

# #NG-14833_Reg
# ${TransDecoder}TransDecoder.LongOrfs -t ${work_dir2}/Trinity.fasta 
# ${TransDecoder}TransDecoder.Predict -t ${work_dir2}/Trinity.fasta 

# #NG-14833_Pas
# ${TransDecoder}TransDecoder.LongOrfs -t ${work_dir3}/Trinity.fasta 
# ${TransDecoder}TransDecoder.Predict -t ${work_dir3}/Trinity.fasta 


###########################################
#Running BLAT
#blat   --minIdentity=80    Trinity.fasta    results_cond1_cond2_k41_coherents_type_0.fa    out_blat_SNP_on_trinity_ID_80.psl

# /bioinfo/BLAT/blat -minIdentity=80 ${work_dir1}/Trinity.fasta ${out_dir1}results_NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1_NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2_Merge_1_Merge_2_NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_k41_coherents_type_0a.fa ${out_dir1}out_blat_SNP_on_trinity_ID_80_Tit.psl
# /bioinfo/BLAT/blat -minIdentity=80 ${work_dir2}/Trinity.fasta ${out_dir2}results_NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1_NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2_Merge_1_Merge_2_Merge_1_Merge_2_Merge_1_Merge_2_NG-14833_Reg_25dpa_2_lib2337k41_coherents_type_0a.fa ${out_dir2}out_blat_SNP_on_trinity_ID_80_Reg.psl
# /bioinfo/BLAT/blat -minIdentity=80 ${work_dir3}/Trinity.fasta ${out_dir3}results_NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1_NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2_NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1_NG-14833_Pas_E_2_lib233712k41_coherents_type_0a.fa ${out_dir3}out_blat_SNP_on_trinity_ID_80_Pas.psl

###########################################
#Running KissDE (optional)

# #!/usr/bin/Rscript
# library(kissDE)

# #NG-14833_Tit
# snp<-kissplice2counts(paste(dir_inout,"results_NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1_NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2_Merge_1_Merge_2_NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_k41_coherents_type_0a.fa",sep=""), pairedEnd=TRUE )
# human_conditions<-c("15dpa","15dpa","15dpa","25dpa","25dpa","25dpa")
# res<-diffExpressedVariants(snp, human_conditions, pvalue=1)
# write.table(res$finalTable, file=paste(dir_inout,"kissDE_Tit_output",sep=""), sep="\t", quote=FALSE)


# #NG-14833_Reg
# snp<-kissplice2counts(paste(dir_inout,"results_NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1_NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2_Merge_1_Merge_2_Merge_1_Merge_2_Merge_1_Merge_2_NG-14833_Reg_25dpa_2_lib2337k41_coherents_type_0a.fa",sep=""), pairedEnd=TRUE )
# human_conditions<-c("15dpa","15dpa","15dpa","25dpa","25dpa","25dpa")
# res<-diffExpressedVariants(snp, human_conditions, pvalue=1)
# write.table(res$finalTable, file=paste(dir_inout,"kissDE_Reg_output",sep=""), sep="\t", quote=FALSE)


# #NG-14833_Pas
# snp<-kissplice2counts(paste(dir_inout,"results_NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1_NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2_NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1_NG-14833_Pas_E_2_lib233712k41_coherents_type_0a.fa",sep=""), pairedEnd=TRUE )
# human_conditions<-c("15dpa","15dpa","15dpa","25dpa","25dpa","25dpa")
# res<-diffExpressedVariants(snp, human_conditions, pvalue=1)
# write.table(res$finalTable, file=paste(dir_inout,"kissDE_Pas_output",sep=""), sep="\t", quote=FALSE)

###########################################
#Running KisSplice2RefTranscriptome

# kissplice2reftranscriptome
# -b Trinity.fasta.transdecoder.bed
# -k result_coherent_type0.fa
# -t out_blat_SNP_on_trinity_ID_80.psl
# -s kissDE_output
# -o mainOutput.tsv
# -l lowQueryCoverageOutput.tsv

# NG-14833_Tit

${kissplice2reftrans} \
-b ${out_dir1_TD}/Trinity.fasta.transdecoder.bed \
-k ${out_dir1}results_NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1_NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2_Merge_1_Merge_2_NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_k41_coherents_type_0a.fa \
-t ${out_dir1}out_blat_SNP_on_trinity_ID_80_Tit.psl \
-s ${out_dir1}kissDE_Tit_output \
-l ${out_dir1}lowQueryCoverageOutput_Tit_mergecodon_pv0.01.tsv \
-o ${out_dir1}kissplice2reftrans_mainOutput_Tit_mergecodon_pv0.01.tsv \
-m ${out_dir1}Merged_bubbles_mergecodon_pv0.01.tsv --merge_codon -pval 0.01

# NG-14833_Reg

${kissplice2reftrans} \
-b ${out_dir2_TD}/Trinity.fasta.transdecoder.bed \
-k ${out_dir2}results_NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1_NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2_Merge_1_Merge_2_Merge_1_Merge_2_Merge_1_Merge_2_NG-14833_Reg_25dpa_2_lib2337k41_coherents_type_0a.fa \
-t ${out_dir2}out_blat_SNP_on_trinity_ID_80_Reg.psl \
-s ${out_dir2}kissDE_Reg_output \
-l ${out_dir2}lowQueryCoverageOutput_Reg_mergecodon_pv0.01.tsv \
-o ${out_dir2}kissplice2reftrans_mainOutput_Reg_mergecodon_pv0.01.tsv \
-m ${out_dir2}Merged_bubbles_mergecodon_pv0.01.tsv --merge_codon -pval 0.01

# NG-14833_Pas

${kissplice2reftrans} \
-b ${out_dir3_TD}/Trinity.fasta.transdecoder.bed \
-k ${out_dir3}results_NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1_NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2_NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1_NG-14833_Pas_E_2_lib233712k41_coherents_type_0a.fa \
-t ${out_dir3}out_blat_SNP_on_trinity_ID_80_Pas.psl \
-s ${out_dir3}kissDE_Pas_output \
-l ${out_dir3}lowQueryCoverageOutput_Pas_mergecodon_pv0.01.tsv \
-o ${out_dir3}kissplice2reftrans_mainOutput_Pas_mergecodon_pv0.01.tsv \
-m ${out_dir3}Merged_bubbles_mergecodon_pv0.01.tsv --merge_codon -pval 0.01


# #########################################

##########
Trinity_path="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5"
work_dir1=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Pas
work_dir2=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Reg
work_dir3=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Tit

# #########################################
# #SuperTranscripts
# #https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts
# #Generate Trinity SuperTranscripts like so:

cd ${work_dir1}
$Trinity_path/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
--trinity_fasta ${work_dir1}/Trinity_old.fasta --incl_malign

cd ${work_dir2}
$Trinity_path/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
--trinity_fasta ${work_dir2}/Trinity_old.fasta --incl_malign

cd ${work_dir3}
$Trinity_path/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
--trinity_fasta ${work_dir3}/Trinity_old.fasta --incl_malign

#######

# #########################################


###############################
# Genes matching to the whole genome of Quinoa
###

#awk 'FNR==NR{a[$2]=$3;next}{print $0,a[$2]?a[$2]:"NA"}' temp2 temp1

id=Tit
file1=/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_${id}/kissplice2reftrans_mainOutput_${id}_mergecodon_pv0.01.tsv
file2=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_${id}/Trinity_old_blastn-max1-ASM168347v1.4_cds.tsv
file3=/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_${id}/kissplice2reftrans_mainOutput_${id}_mergecodon_pv0.01_blastn_ASM168347v1_cds.tsv

awk 'FNR==NR{a[$1]=$2;next}{print $0,"\t",a[$1]?a[$1]:"NA"}' $file2 $file1 > $file3


###############################
# Genes matching to the whole genome of Arobidopsis 
###

file1=/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_${id}/kissplice2reftrans_mainOutput_${id}_mergecodon_pv0.01_blastn_ASM168347v1_cds.tsv
file2=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_${id}/Trinity_old_blastn-max1-000001735.4_cds.tsv
file3=/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_${id}/kissplice2reftrans_mainOutput_${id}_mergecodon_pv0.01_blastn_ASM168347v1_000001735.4_cds.tsv

awk 'FNR==NR{a[$1]=$2;next}{print $0,"\t",a[$1]?a[$1]:"NA"}' $file2 $file1 > $file3



echo "Script done...."

