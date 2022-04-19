#!/bin/bash


#Rna-seq using STAR
#./2.1.2_Trinity_Genome-guided_Master-ass.sh > 2.1.2_Trinity_Genome-guided_Master-ass.sh.log 2>&1


echo "run script for rna-seq-analysis"
##########################################
##########################################
INDIR=/home/gala0002/proj/proj_dir/NG-14833_1.1_sort-trim/
REF_Quinoa=/home/gala0002/proj/proj_dir/REF_Genome/Ref_Chenopodium_Quinoa/
work_dir_GG_master=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_GG_Master_Tit-Reg-Pas/
Trinity="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5/Trinity"
Trinity_path="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5"
proj_dir="/home/gala0002/proj/proj_dir/"

GG_Script=$Trinity_path/Analysis/SuperTranscripts/AllelicVariants/2.1.2_Trinity_Genome-guided.py

#Genome Guded De novo assembly, using Trinity

#######
ulimit -n 65535

# Master GG de-novo assembly
#####################
#Tit-Reg-Pas_15dpa-25dpa
#####################

# nice -n 5 $GG_Script \
# --st_fa ${REF_Quinoa}GCF_001683475.1_ASM168347v1_genomic.fna \
# --st_gtf ${REF_Quinoa}GCF_001683475.1_ASM168347v1_genomic.gff \
# -p ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_2.fq.gz \
# -o ${work_dir_GG_master}GG-Master-ass_Tit-Reg-Pas_15dpa-25dpa -m 300000000000 -t 35


# ##################### didn't do this #####

# #Tit-15dpa and #Tit-25dpa
# nice -n 5 $GG_Script \
# --st_fa ${REF_Quinoa}GCF_001683475.1_ASM168347v1_genomic.fna \
# --st_gtf ${REF_Quinoa}GCF_001683475.1_ASM168347v1_genomic.gff \
# -p ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_1.fq.gz ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_2.fq.gz \
# -o ${work_dir_GG_master}GG-Master-ass_Tit_15dpa-25dpa -m 199000000000 -t 35

# #######

# #Reg-15dpa and #Reg-25dpa
# nice -n 5 $GG_Script \
# --st_fa ${REF_Quinoa}GCF_001683475.1_ASM168347v1_genomic.fna \
# --st_gtf ${REF_Quinoa}GCF_001683475.1_ASM168347v1_genomic.gff \
# -p ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_2.fq.gz \
# -o ${work_dir_GG_master}GG-Master-ass_Reg_15dpa-25dpa -m 199000000000 -t 35

# #######

# #Pas-15dpa and #Pas-25dpa
# nice -n 5 $GG_Script \
# --st_fa ${REF_Quinoa}GCF_001683475.1_ASM168347v1_genomic.fna \
# --st_gtf ${REF_Quinoa}GCF_001683475.1_ASM168347v1_genomic.gff \
# -p ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_2.fq.gz \
# -o ${work_dir_GG_master}GG-Master-ass_Pas_15dpa-25dpa -m 199000000000 -t 35

# ##################### didn't do this #####

##########################################
#####################
##########################################
#2. Generating a Trinity Genome guided de-novo RNA-Seq assembly

# nice -n 5 $Trinity --genome_guided_bam ${work_dir_GG_master}GG-Master-ass_Tit-Reg-Pas_15dpa-25dpa/Aligned.sortedByCoord.out.bam \
# --genome_guided_max_intron 10000 --SS_lib_type RF \
# --max_memory 200G --CPU 30


# #########################################
# ### Trinity Stats
# #2.1. Evaluating the quality of the assembly

# cd ${work_dir_GG_master}
# $Trinity_path/util/TrinityStats.pl ${work_dir_GG_master}trinity_out_dir/Trinity-GG.fasta >& ${work_dir_GG_master}stats

# # #########################################
# # #SuperTranscripts
# # #https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts
# # #Generate Trinity SuperTranscripts like so:

# $Trinity_path/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
#  --trinity_fasta ${work_dir_GG_master}trinity_out_dir/Trinity-GG.fasta --incl_malign

# #########################################
# #Variant Calling
# https://github.com/trinityrnaseq/trinityrnaseq/wiki/Variant-Calling


# #########################################
# #Variant Calling
# https://github.com/trinityrnaseq/trinityrnaseq/wiki/Variant-Calling
##########################################
#####################

#######
ulimit -n 65535
ulimit -c unlimited

#Tit-15dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_GG_master}trinity_genes.fasta \
--st_gtf ${work_dir_GG_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_2.fq.gz \
-o ${work_dir_GG_master}variants_Tit-15dpa -m 250000000000 -t 45 > ${work_dir_GG_master}variants_Tit-15dpa.log 2>&1


#Tit-25dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_GG_master}trinity_genes.fasta \
--st_gtf ${work_dir_GG_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_1.fq.gz ${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_2.fq.gz \
-o ${work_dir_GG_master}variants_Tit-25dpa -m 250000000000 -t 45 > ${work_dir_GG_master}variants_Tit-25dpa.log 2>&1


#######
ulimit -n 65535
ulimit -c unlimited

#Reg-15dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_GG_master}trinity_genes.fasta \
--st_gtf ${work_dir_GG_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_1.fq.gz ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_2.fq.gz \
-o ${work_dir_GG_master}variants_Reg-15dpa -m 250000000000 -t 45 > ${work_dir_GG_master}variants_Reg-15dpa.log 2>&1


#Reg-25dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_GG_master}trinity_genes.fasta \
--st_gtf ${work_dir_GG_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_2.fq.gz \
-o ${work_dir_GG_master}variants_Reg-25dpa -m 250000000000 -t 45 > ${work_dir_GG_master}variants_Reg-25dpa.log 2>&1


#######
ulimit -n 65535
ulimit -c unlimited

#Pas-15dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_GG_master}trinity_genes.fasta \
--st_gtf ${work_dir_GG_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_2.fq.gz \
-o ${work_dir_GG_master}variants_Pas-15dpa -m 250000000000 -t 45 > ${work_dir_GG_master}variants_Pas-15dpa.log 2>&1


#Pas-25dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_GG_master}trinity_genes.fasta \
--st_gtf ${work_dir_GG_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_2.fq.gz \
-o ${work_dir_GG_master}variants_Pas-25dpa -m 250000000000 -t 45 > ${work_dir_GG_master}variants_Pas-25dpa.log 2>&1


#######



############################################
# 1.0 Differential Transcript Usage via SuperTranscripts 
# DiffTranscriptUsage (DTU)
############################################
#https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage

$Trinity_path/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
--genes_fasta ${work_dir_GG_master}trinity_genes.fasta \
--genes_gtf ${work_dir_GG_master}trinity_genes.gtf \
--samples_file $proj_dir/Samples_Master_Tit-Reg-Pas.txt \
--out_prefix DTU --aligner STAR > SuperTranscripts-DTU.log 2>&1


############################################
# 2.1. Evaluating the quality of the assembly
# Assembly metrics
############################################
#https://southgreenplatform.github.io/tutorials//bioanalysis/trinity/

perl $Trinity_path/util/align_and_estimate_abundance.pl \
--transcripts Trinity-GG.fasta \
--seqType fq \
--samples_file $proj_dir/Samples_Master_Tit-Reg-Pas.txt \
--est_method salmon \
--trinity_mode --prep_reference \
--output_dir outdir_estimate-ab > estimate-ab.log 2>&1

############################################
# genes 

$Trinity_path/util/abundance_estimates_to_matrix.pl --est_method salmon \
--gene_trans_map none \
--out_prefix Trinity_genes \
--name_sample_by_basedir \
Tit_15dpa_rep1/quant.sf.genes \
Tit_15dpa_rep2/quant.sf.genes \
Tit_15dpa_rep3/quant.sf.genes \
Tit_25dpa_rep1/quant.sf.genes \
Tit_25dpa_rep2/quant.sf.genes \
Tit_25dpa_rep3/quant.sf.genes \
Reg_15dpa_rep1/quant.sf.genes \
Reg_15dpa_rep2/quant.sf.genes \
Reg_15dpa_rep3/quant.sf.genes \
Reg_25dpa_rep1/quant.sf.genes \
Reg_25dpa_rep2/quant.sf.genes \
Reg_25dpa_rep3/quant.sf.genes \
Pas_15dpa_rep1/quant.sf.genes \
Pas_15dpa_rep2/quant.sf.genes \
Pas_15dpa_rep3/quant.sf.genes \
Pas_25dpa_rep1/quant.sf.genes \
Pas_25dpa_rep2/quant.sf.genes \
Pas_25dpa_rep3/quant.sf.genes > Trinity_genes_salmon.log 2>&1

#trans 
$Trinity_path/util/abundance_estimates_to_matrix.pl --est_method salmon \
--gene_trans_map none \
--out_prefix Trinity_trans \
--name_sample_by_basedir \
Tit_15dpa_rep1/quant.sf \
Tit_15dpa_rep2/quant.sf \
Tit_15dpa_rep3/quant.sf \
Tit_25dpa_rep1/quant.sf \
Tit_25dpa_rep2/quant.sf \
Tit_25dpa_rep3/quant.sf \
Reg_15dpa_rep1/quant.sf \
Reg_15dpa_rep2/quant.sf \
Reg_15dpa_rep3/quant.sf \
Reg_25dpa_rep1/quant.sf \
Reg_25dpa_rep2/quant.sf \
Reg_25dpa_rep3/quant.sf \
Pas_15dpa_rep1/quant.sf \
Pas_15dpa_rep2/quant.sf \
Pas_15dpa_rep3/quant.sf \
Pas_25dpa_rep1/quant.sf \
Pas_25dpa_rep2/quant.sf \
Pas_25dpa_rep3/quant.sf > Trinity_trans_salmon.log 2>&1


#Compute N50 based on the top-most highly expressed transcripts (Ex50)
$Trinity_path/util/misc/contig_ExN50_statistic.pl Trinity_trans.isoform.TMM.EXPR.matrix Trinity-GG.fasta > ExN50.stats

############################################
# 2.2. Identifying differentially expressed (DE) transcripts
############################################
#Extracting differentially expressed transcripts and generating heatmaps

############################################
$Trinity_path/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix Trinity_trans.isoform.counts.matrix \
--samples_file $proj_dir/Samples_Master_Tit-Reg-Pas.txt \
--method DESeq2 \
--output DESeq2_trans > Trinity_DE-run-trans.log 2>&1

#Extracting differentially expressed transcripts and generating heatmaps
#Extract those differentially expressed (DE) transcripts that are at least 4-fold (C is set to 2^(2) ) differentially expressed at a significance of <= 0.001 (-P 1e-3) in any of the pairwise sample comparisons

cd DESeq2_trans/
$Trinity_path/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../Trinity_trans.isoform.TMM.EXPR.matrix \
--samples $proj_dir/Samples_Master_Tit-Reg-Pas.txt -P 1e-3 -C 2 > Trinity_DE-analyze-trans.log 2>&1
cd ..

############################################
#Run the DE analysis at the gene level

$Trinity_path/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix Trinity_genes.isoform.counts.matrix \
--samples_file $proj_dir/Samples_Master_Tit-Reg-Pas.txt \
--method DESeq2 \
--output DESeq2_gene > Trinity_DGE-run-genes.log 2>&1

#Extracting differentially expressed transcripts and generating heatmaps
#Extract those differentially expressed (DE) transcripts that are at least 4-fold (C is set to 2^(2) ) differentially expressed at a significance of <= 0.001 (-P 1e-3) in any of the pairwise sample comparisons
cd DESeq2_gene/
$Trinity_path/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../Trinity_genes.isoform.TMM.EXPR.matrix \
--samples $proj_dir/Samples_Master_Tit-Reg-Pas.txt -P 1e-3 -C 2 > Trinity_DGE-analyze-genes.log 2>&1
cd ..

############################################

######################### ---OPTIONAL ----

$Trinity_path/Analysis/SuperTranscripts/AllelicVariants/VCF_to_annotated_SNP_report.pl \
--gff3 trinity_genes.gtf \
--genome trinity_genes.fasta \
--vcf /home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_GG_Master_Tit-Reg-Pas/variants_Pas-15dpa/filtered_output.vcf > VCF_to_annotated_SNP_report.log 2>&1
######################### ---OPTIONAL ----

############################################




echo "Script done...."

#######



