#!/bin/bash


#Rna-seq using STAR
#./4.0_DET.sh > 4.0_DET.sh.log 2>&1


echo "run script for rna-seq-analysis"

############################################
# 2.2. Identifying differentially expressed (DE) transcripts
############################################
#Extracting differentially expressed transcripts and generating heatmaps

Trinity="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5/Trinity"
Trinity_path="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5"
work_dir=/home/gala0002/proj/proj_dir/
#ref=/home/gala0002/proj/proj_Ramesh/Ref_Nicotiana_benthamiana/

cd $work_dir


############################################
#Build Transcript and Gene Expression Matrices
#https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#express-output
# perl $Trinity_path/util/align_and_estimate_abundance.pl \
# --transcripts Trinity.fasta \
# --seqType fq \
# --samples_file Metadata_CDS_Mock-PMTVWT-delta8k_merge.txt \
# --est_method salmon \
# --trinity_mode --prep_reference \
# --output_dir outdir_estimate-ab > estimate-ab.log 2>&1

############################################
#Run the DE analysis at the gene level

#DESeq2
#for i in Pas-Reg-Tit Pas-vs-Reg-Tit
#Pas-Reg-Tit Pas-vs-Reg-Tit Reg-vs-Pas-Tit Tit-vs-Pas-Reg
#Pas-vs-Reg Pas-vs-Tit

for i in Pas-vs-Reg
do 

echo ${i}

mkdir -p ${work_dir}"DESeq2_genes_Quinoa_"${i}

cd $work_dir

<<COMMENT
COMMENT

$Trinity_path/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix Express_counts_allsamp.matrix \
--samples_file Metadata_${i}.txt \
--method DESeq2 \
--output DESeq2_genes_Quinoa_${i} > DGE-run_Quinoa_${i}.log 2>&1


#Extracting differentially expressed transcripts and generating heatmaps
#Extract those differentially expressed (DE) transcripts that are at least 4-fold (C is set to 2^(2) ) differentially expressed at a significance of <= 0.001 (-P 1e-3) in any of the pairwise sample comparisons
#-C 1.0

cd DESeq2_genes_Quinoa_${i}/
nice -n 5 $Trinity_path/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../Express_counts_allsamp.matrix \
--samples ../Metadata_${i}.txt -P 1e-2 -C 1 > DGE-analyze_Quinoa_${i}.log 2>&1

rename 's/Express_counts_allsamp.matrix.//' *

#rm *_vs_P*

#rm diffExpr.P5e-2_C0.matrix.RData

done

####

############################################

echo "Script done...."


