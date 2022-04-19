#!/bin/bash -l


echo "run script for rna-seq-analysis Cleaning and QC"

#./1.1_rRNA-and-AdaptRM.sh > 1.1_rRNA-and-AdaptRM.sh.log 2>&1

#1. Raw Data rRNA and AdaptRM Assessment

work_dir=/home/proj/

mkdir -p ${work_dir}NG-14833_1.1_sort-trim/
out_dir1=${work_dir}NG-14833_1.1_sort-trim/

################################
#####---Lane: NG-14833
################################

cd ${work_dir}Asa-Grimberg_Quinoa-GATC-2018/
#for nbr in Sample_484-10-1/; do
#for nbr in NG-14833_*/; do
#nbr=$(echo $nbr | sed 's=/[^/]*$==;s/\.$//')
nbr=$1
mkdir ${out_dir1}${nbr}

#----
SORTMERNADIR=/bioinfo/sortmerna-2.1b
TRIMMOMATIC=/bioinfo/Trimmomatic-0.36

data_dir=${work_dir}Asa-Grimberg_Quinoa-GATC-2018/
scripts=/bioinfo/sortmerna-2.1b/scripts/

temp_dir=${out_dir1}${nbr}"/"
cd ${temp_dir}

echo ${temp_dir}, ${out_dir1}${nbr}, $data_dir, $nbr

gunzip -c ${data_dir}${nbr}"_1.fastq.gz" > ${temp_dir}${nbr}_R1.fq
#gunzip -c ${data_dir}${nbr}"_2.fastq.gz" > ${temp_dir}${nbr}_R2.fq

