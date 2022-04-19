#!/bin/bash -l

echo "run script for rna-seq-analysis Cleaning and QC"

#./1.2_Merge-lanes.sh > 1.2_Merge-lanes.sh.log 2>&1

#1. Raw Data QC Assessment
work_dir=/home/proj/NG-14833_1.1_sort-trim/
cd ${work_dir}

##############--NG-14833_Reg_15dpa_2
#Command to Merge two lanes
mkdir -p ${work_dir}NG-14833_Reg_15dpa_2_lib233718_Merge/
out_dir1=${work_dir}NG-14833_Reg_15dpa_2_lib233718_Merge/

SAMP1="NG-14833_Reg_15dpa_2_lib233718_5795_1"
SAMP2="NG-14833_Reg_15dpa_2_lib233718_5767_1"

nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L1_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L1_2.fq

nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L2_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L2_2.fq

cat ${out_dir1}L1_1.fq ${out_dir1}L2_1.fq > ${out_dir1}Merge_1.fq
cat ${out_dir1}L1_2.fq ${out_dir1}L2_2.fq > ${out_dir1}Merge_2.fq

#Command to run fastqc
nice -n 5 fastqc -o ${out_dir1} -t 15 --noextract ${out_dir1}Merge_1.fq ${out_dir1}Merge_2.fq

#remove temp files
rm ${out_dir1}L*_*.fq
#rm ${out_dir1}Merge_*.fq
##############

##############--NG-14833_Reg_15dpa_3
#Command to Merge two lanes
mkdir -p ${work_dir}NG-14833_Reg_15dpa_3_lib233719_Merge/
out_dir1=${work_dir}NG-14833_Reg_15dpa_3_lib233719_Merge/

SAMP1="NG-14833_Reg_15dpa_3_lib233719_5767_1"
SAMP2="NG-14833_Reg_15dpa_3_lib233719_5795_1"

nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L1_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L1_2.fq

nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L2_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L2_2.fq

cat ${out_dir1}L1_1.fq ${out_dir1}L2_1.fq > ${out_dir1}Merge_1.fq
cat ${out_dir1}L1_2.fq ${out_dir1}L2_2.fq > ${out_dir1}Merge_2.fq

#Command to run fastqc
nice -n 5 fastqc -o ${out_dir1} -t 15 --noextract ${out_dir1}Merge_1.fq ${out_dir1}Merge_2.fq

#remove temp files
rm ${out_dir1}L*_*.fq
#rm ${out_dir1}Merge_*.fq

##############

##############--NG-14833_Reg_25dpa_1
#Command to Merge two lanes
mkdir -p ${work_dir}NG-14833_Reg_25dpa_1_lib233720_Merge/
out_dir1=${work_dir}NG-14833_Reg_25dpa_1_lib233720_Merge/

SAMP1="NG-14833_Reg_25dpa_1_lib233720_5767_1"
SAMP2="NG-14833_Reg_25dpa_1_lib233720_5795_1"

nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L1_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L1_2.fq

nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L2_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L2_2.fq

cat ${out_dir1}L1_1.fq ${out_dir1}L2_1.fq > ${out_dir1}Merge_1.fq
cat ${out_dir1}L1_2.fq ${out_dir1}L2_2.fq > ${out_dir1}Merge_2.fq

#Command to run fastqc
nice -n 5 fastqc -o ${out_dir1} -t 15 --noextract ${out_dir1}Merge_1.fq ${out_dir1}Merge_2.fq

#remove temp files
rm ${out_dir1}L*_*.fq
#rm ${out_dir1}Merge_*.fq

##############


##############--NG-14833_Tit_15dpa_2_lib233706
#Command to Merge two lanes
mkdir -p ${work_dir}NG-14833_Tit_15dpa_2_lib233706_Merge/
out_dir1=${work_dir}NG-14833_Tit_15dpa_2_lib233706_Merge/

SAMP1="NG-14833_Tit_15dpa_2_lib233706_5747_1"
SAMP2="NG-14833_Tit_15dpa_2_lib233706_5767_7"

nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L1_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L1_2.fq

nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L2_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L2_2.fq

cat ${out_dir1}L1_1.fq ${out_dir1}L2_1.fq > ${out_dir1}Merge_1.fq
cat ${out_dir1}L1_2.fq ${out_dir1}L2_2.fq > ${out_dir1}Merge_2.fq

#Command to run fastqc
nice -n 5 fastqc -o ${out_dir1} -t 15 --noextract ${out_dir1}Merge_1.fq ${out_dir1}Merge_2.fq

#remove temp files
rm ${out_dir1}L*_*.fq
#rm ${out_dir1}Merge_*.fq

##############


##############--NG-14833_Tit_25dpa_3_lib233710
#Command to Merge two lanes
mkdir -p ${work_dir}NG-14833_Tit_25dpa_3_lib233710_Merge/
out_dir1=${work_dir}NG-14833_Tit_25dpa_3_lib233710_Merge/

SAMP1="NG-14833_Tit_25dpa_3_lib233710_5747_1"
SAMP2="NG-14833_Tit_25dpa_3_lib233710_5767_7"

nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L1_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP1}/${SAMP1}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L1_2.fq

nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_1.fq.gz > ${out_dir1}L2_1.fq
nice -n 5 gunzip -c ${work_dir}${SAMP2}/${SAMP2}-sortmerna-trimmomatic_2.fq.gz > ${out_dir1}L2_2.fq

cat ${out_dir1}L1_1.fq ${out_dir1}L2_1.fq > ${out_dir1}Merge_1.fq
cat ${out_dir1}L1_2.fq ${out_dir1}L2_2.fq > ${out_dir1}Merge_2.fq

#Command to run fastqc
nice -n 5 fastqc -o ${out_dir1} -t 15 --noextract ${out_dir1}Merge_1.fq ${out_dir1}Merge_2.fq

#remove temp files
rm ${out_dir1}L*_*.fq
#rm ${out_dir1}Merge_*.fq

##############

echo "Done script..."


