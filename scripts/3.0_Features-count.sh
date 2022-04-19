#!/bin/bash


#Separate reads
##############################################################################################################

#./3.0_Features-count.sh > 3.0_Features-count.sh.log 2>&1

echo "run script for rna-seq-analysis"

######################
work_dir=/home/gala0002/proj/proj_dir/

mkdir -p ${work_dir}NG-14833_3.0_Salmon_Quinoa_Genome/
out_dir=${work_dir}NG-14833_3.0_Salmon_Quinoa_Genome/

in_dir=${work_dir}NG-14833_2.0_Align-STAR_Quinoa_Genome/

cd ${in_dir}

for nbr in `ls $in_dir`
do

nbr1=$(echo $nbr | sed 's=/[^/]*$==;s/\.$//')
mkdir -p ${out_dir}${nbr}/

temp_dir=${out_dir}${nbr}/
cd ${temp_dir}

#All-reads
echo "Processing sample: ${in_dir}${nbr}/${nbr1}"

######## #Salmon #####
#https://salmon.readthedocs.io/en/latest/salmon.html
#conda create -n salmon salmon

conda activate salmon

#Ref_trans="/home/gala0002/proj/proj_Sophie-Laura/REF/Altso1_AssemblyScaffolds/Altso1_AssemblyScaffolds_Repeatmasked_trans.fa"
Ref_trans="/home/gala0002/proj/proj_dir/REF_Genome/Ref_Chenopodium_Quinoa/GCF_001683475.1_ASM168347v1_trans.fa"

echo "Processing sample ${i}"
salmon quant -t ${Ref_trans} -p 60 -l A \
-a ${in_dir}${nbr}/${nbr1}-sort-trim-STARAligned.toTranscriptome.out.bam \
-o ${out_dir}${nbr1}/salmon

conda deactivate

<<COMMENT
#######################################################
COMMENT


done

############

echo "Script done all...."



