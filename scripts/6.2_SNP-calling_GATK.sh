#!/bin/bash


#Rna-seq using STAR
#./6.2_SNP-calling_GATK.sh > 6.2_SNP-calling_GATK.sh.log 2>&1

#https://github.com/Nourolah/Quinoa-RNASeq---Variant-Calling/wiki/Variant-Calls-on-RNASeq

echo "run script for rna-seq-analysis"
##########################################
##########################################
INDIR=/home/gala0002/proj/proj_dir/NG-14833_1.1_sort-trim
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

###########################################
#De novo assembly of reads using Trinity










echo "Script done...."


