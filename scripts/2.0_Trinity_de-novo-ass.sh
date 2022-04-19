#!/bin/bash


#Rna-seq using STAR
#./2.0_Trinity_de-novo-ass.sh > 2.0_Trinity_de-novo-ass.sh.log 2>&1

echo "run script for rna-seq-analysis"
##########################################
##########################################
INDIR=/home/gala0002/proj/proj_dir/NG-14833_1.1_sort-trim/
work_dir1=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Tit/
work_dir2=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Reg/
work_dir3=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Pas/
Trinity="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5/Trinity"
Trinity_path="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5/"
proj_dir="/home/gala0002/proj/proj_dir/"

##########################################
#De novo assembly of reads using Trinity
#${TRINITY_HOME}/Trinity --seqType fq --SS_lib_type RF  \
#--left RNASEQ_data/Sp_log.left.fq.gz,RNASEQ_data/Sp_hs.left.fq.gz,RNASEQ_data/Sp_ds.left.fq.gz,RNASEQ_data/Sp_plat.left.fq.gz \
#--right RNASEQ_data/Sp_log.right.fq.gz,RNASEQ_data/Sp_hs.right.fq.gz,RNASEQ_data/Sp_ds.right.fq.gz,RNASEQ_data/Sp_plat.right.fq.gz \
#--CPU 2 --max_memory 1G

#MergeDIR1=/home/gala0002/proj/proj_dir/NG-14833_1.1_sort-trim/NG-14833_Tit_15dpa_2_lib233706_Merge/
#MergeDIR2=/home/gala0002/proj/proj_dir/NG-14833_1.1_sort-trim/NG-14833_Tit_25dpa_3_lib233710_Merge/

#gzip -c Merge_1.fq > Merge_1.fq.gz
#gzip -c ${MergeDIR1}Merge_2.fq > ${MergeDIR1}Merge_2.fq.gz
#gzip -c ${MergeDIR2}Merge_1.fq > ${MergeDIR2}Merge_1.fq.gz
#gzip -c ${MergeDIR2}Merge_2.fq > ${MergeDIR2}Merge_2.fq.gz

##########################################
#1.4. Normalisation using Trinity
#ref https://southgreenplatform.github.io/tutorials//bioanalysis/trinity/

# maindir="/home/gala0002/proj/proj_dir/"
# mkdir -p ${maindir}"NG-14833_2.0_Trinity_Tit_norm"

# perl /data/bioinfo/trinityrnaseq-Trinity-v2.6.5/util/insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
# --left ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_1.fq.gz \
# --right ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_2.fq.gz \
# --pairs_together --PARALLEL_STATS --CPU 25 --output ${maindir}"NG-14833_2.0_Trinity_Tit_norm"/

##########################################
#2. Generating a Trinity de novo RNA-Seq assembly

#--left ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_1.fq.gz \
#--right ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_2.fq.gz \

#Trinity --seqType fq --max_memory 80G --CPU 8 --normalize_by_read_set --samples_file samples.txt --output trinity_OUT
# nice -n 5 ${Trinity} --seqType fq --SS_lib_type RF \
# --CPU 25 --max_memory 80G --normalize_by_read_set --samples_file Samples_Tit.txt --output ${work_dir} --workdir ${work_dir}

nice -n 5 ${Trinity} --seqType fq --SS_lib_type RF \
--CPU 25 --max_memory 80G --normalize_by_read_set --samples_file ${proj_dir}Samples_Reg.txt --output ${work_dir2} --workdir ${work_dir2}

nice -n 5 ${Trinity} --seqType fq --SS_lib_type RF \
--CPU 25 --max_memory 80G --normalize_by_read_set --samples_file ${proj_dir}Samples_Pas.txt --output ${work_dir3} --workdir ${work_dir3}


# #########################################
# ### Trinity Stats
# #2.1. Evaluating the quality of the assembly
# cd ${work_dir1}
# $Trinity_path/util/TrinityStats.pl ${work_dir1}Trinity.fasta >& stats
cd ${work_dir2}
$Trinity_path/util/TrinityStats.pl ${work_dir2}Trinity.fasta >& stats
cd ${work_dir3}
$Trinity_path/util/TrinityStats.pl ${work_dir3}Trinity.fasta >& stats

# #########################################
# #Reads mapping back rate :
# #A typical ‘good’ assembly has ~80 % reads mapping to the assembly and ~80% are properly paired
# #Alignment methods : bowtie2 -RSEM, kallisto, salmon --est_method

#mkdir -p ${work_dir1}"out_align_and_estimate_abundance"
mkdir -p ${work_dir2}"out_align_and_estimate_abundance"
mkdir -p ${work_dir3}"out_align_and_estimate_abundance"

#outdir1=${work_dir1}"out_align_and_estimate_abundance"
outdir2=${work_dir2}"out_align_and_estimate_abundance"
outdir3=${work_dir3}"out_align_and_estimate_abundance"

# perl $Trinity_path/util/align_and_estimate_abundance.pl \
# --transcripts ${work_dir1}Trinity.fasta \
# --seqType fq \
# --samples_file ${proj_dir}Samples_Tit.txt \
# --est_method salmon \
# --trinity_mode --prep_reference \
# --output_dir $outdir1

perl $Trinity_path/util/align_and_estimate_abundance.pl \
--transcripts ${work_dir2}Trinity.fasta \
--seqType fq \
--samples_file ${proj_dir}Samples_Reg.txt \
--est_method salmon \
--trinity_mode --prep_reference \
--output_dir $outdir2

perl $Trinity_path/util/align_and_estimate_abundance.pl \
--transcripts ${work_dir3}Trinity.fasta \
--seqType fq \
--samples_file ${proj_dir}Samples_Pas.txt \
--est_method salmon \
--trinity_mode --prep_reference \
--output_dir $outdir3 

# #########################################

echo "Script done...."
