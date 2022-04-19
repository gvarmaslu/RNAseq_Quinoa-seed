#!/bin/bash


#Rna-seq using STAR
#./2.0_Trinity_de-novo_Master-ass.sh > 2.0_Trinity_de-novo_Master-ass.sh.log 2>&1

echo "run script for rna-seq-analysis"
##########################################
##########################################
INDIR=/home/gala0002/proj/proj_dir/NG-14833_1.1_sort-trim/
work_dir_master=/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Master_Tit-Reg-Pas/
#mkdir -p $work_dir_master
Trinity="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5/Trinity"
Trinity_path="/data/bioinfo/trinityrnaseq-Trinity-v2.6.5"
proj_dir="/home/gala0002/proj/proj_dir/"

##########################################
#De novo assembly of reads using Trinity
#${TRINITY_HOME}/Trinity --seqType fq --SS_lib_type RF  \
#--left RNASEQ_data/Sp_log.left.fq.gz,RNASEQ_data/Sp_hs.left.fq.gz,RNASEQ_data/Sp_ds.left.fq.gz,RNASEQ_data/Sp_plat.left.fq.gz \
#--right RNASEQ_data/Sp_log.right.fq.gz,RNASEQ_data/Sp_hs.right.fq.gz,RNASEQ_data/Sp_ds.right.fq.gz,RNASEQ_data/Sp_plat.right.fq.gz \
#--CPU 2 --max_memory 1G

##########################################
#2. Generating a Trinity de novo RNA-Seq assembly

# nice -n 5 ${Trinity} --seqType fq --SS_lib_type RF \
# --CPU 25 --max_memory 100G --normalize_by_read_set --samples_file ${proj_dir}Samples_Master_Tit-Reg-Pas.txt --output ${work_dir_master} --workdir ${work_dir_master}

# #########################################
# ### Trinity Stats
# #2.1. Evaluating the quality of the assembly
# cd ${work_dir1}
# $Trinity_path/util/TrinityStats.pl ${work_dir1}Trinity.fasta >& stats

# cd ${work_dir_master}
# $Trinity_path/util/TrinityStats.pl ${work_dir_master}Trinity.fasta >& ${work_dir_master}stats

# #########################################
# #SuperTranscripts
# #https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts
# #Generate Trinity SuperTranscripts like so:

# $Trinity_path/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
# --trinity_fasta ${work_dir_master}Trinity.fasta --incl_malign

# #########################################
# #Variant Calling
# https://github.com/trinityrnaseq/trinityrnaseq/wiki/Variant-Calling

#######
ulimit -n 65535

#Tit-15dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_master}trinity_genes.fasta \
--st_gtf ${work_dir_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Tit_15dpa_1_lib233705_5747_1/NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_2_lib233706_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Tit_15dpa_3_lib233707_5747_1/NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_2.fq.gz \
-o ${work_dir_master}variants_Tit-15dpa -m 199000000000 -t 45


#Tit-25dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_master}trinity_genes.fasta \
--st_gtf ${work_dir_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_1.fq.gz ${INDIR}NG-14833_Tit_25dpa_1_lib233708_5747_1/NG-14833_Tit_25dpa_1_lib233708_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_2_lib233709_5747_1/NG-14833_Tit_25dpa_2_lib233709_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Tit_25dpa_3_lib233710_Merge/Merge_2.fq.gz \
-o ${work_dir_master}variants_Tit-25dpa -m 199000000000 -t 45


#######

#Reg-15dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_master}trinity_genes.fasta \
--st_gtf ${work_dir_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_1.fq.gz ${INDIR}NG-14833_Reg_15dpa_1_lib233717_5767_1/NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_2_lib233718_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_15dpa_3_lib233719_Merge/Merge_2.fq.gz \
-o ${work_dir_master}variants_Reg-15dpa -m 199000000000 -t 45


#Reg-25dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_master}trinity_genes.fasta \
--st_gtf ${work_dir_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Reg_25dpa_1_lib233720_Merge/Merge_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_2_lib233721_5747_6/NG-14833_Reg_25dpa_2_lib233721_5747_6-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Reg_25dpa_3_lib233722_5747_6/NG-14833_Reg_25dpa_3_lib233722_5747_6-sortmerna-trimmomatic_2.fq.gz \
-o ${work_dir_master}variants_Reg-25dpa -m 199000000000 -t 45



#######

#Pas-15dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_master}trinity_genes.fasta \
--st_gtf ${work_dir_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Pas_E_1_lib233711_5747_1/NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_2_lib233712_5747_1/NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_E_3_lib233713_5767_1/NG-14833_Pas_E_3_lib233713_5767_1-sortmerna-trimmomatic_2.fq.gz \
-o ${work_dir_master}variants_Pas-15dpa -m 199000000000 -t 45


#Pas-25dpa
nice -n 5 $Trinity_path/Analysis/SuperTranscripts/AllelicVariants/run_variant_calling.py \
--st_fa ${work_dir_master}trinity_genes.fasta \
--st_gtf ${work_dir_master}trinity_genes.gtf \
-p ${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_1.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_1.fq.gz ${INDIR}NG-14833_Pas_M_1_lib233714_5767_1/NG-14833_Pas_M_1_lib233714_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_2_lib233715_5767_1/NG-14833_Pas_M_2_lib233715_5767_1-sortmerna-trimmomatic_2.fq.gz,${INDIR}NG-14833_Pas_M_3_lib233716_5767_1/NG-14833_Pas_M_3_lib233716_5767_1-sortmerna-trimmomatic_2.fq.gz \
-o ${work_dir_master}variants_Pas-25dpa -m 199000000000 -t 45




# #########################################
# #Reads mapping back rate :
# #A typical ‘good’ assembly has ~80 % reads mapping to the assembly and ~80% are properly paired
# #Alignment methods : bowtie2 -RSEM, kallisto, salmon --est_method
# #########################################

ulimit -n 65535
ulimit -c unlimited

/usr/lib/jvm/java-8-openjdk-amd64/bin/java -Xmx100g -jar /bioinfo/GATK/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar \
-T SplitNCigarReads -R /data/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Master_Tit-Reg-Pas/trinity_genes.fasta \
-I dedupped.valid.bam -o splitNCigar.bam  -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS --validation_strictness LENIENT > GATKrun-splitN.log 2>&1

touch trinity_genes.fasta.gatk_chkpts/splitNCigarReads.ok

/usr/lib/jvm/java-8-openjdk-amd64/bin/java -Xmx100g -jar /bioinfo/GATK/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar \
-T HaplotypeCaller -R /data/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Master_Tit-Reg-Pas/trinity_genes.fasta \
-I ./splitNCigar.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.vcf > GATKrun-Hapcall.log 2>&1

touch trinity_genes.fasta.gatk_chkpts/haplotypecaller.ok 

/usr/lib/jvm/java-8-openjdk-amd64/bin/java -Xmx100g -jar /bioinfo/GATK/GenomeAnalysisTK-3.8-1/GenomeAnalysisTK.jar \
-T VariantFiltration -R /data/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Master_Tit-Reg-Pas/trinity_genes.fasta \
-V ./output.vcf -window 35 -cluster 3 \
--filterName "FS" -filter "FS > 30.0" \
--filterName "QD" -filter "QD < 2.0" -o filtered_output.vcf > GATKrun-VarFil.log 2>&1

touch trinity_genes.fasta.gatk_chkpts/variant_filt.ok



############################################
# 1.0 Differential Transcript Usage via SuperTranscripts 
# DiffTranscriptUsage (DTU)
############################################
#https://github.com/trinityrnaseq/trinityrnaseq/wiki/DiffTranscriptUsage

$Trinity_path/Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl \
--genes_fasta ${work_dir_master}trinity_genes.fasta \
--genes_gtf ${work_dir_master}trinity_genes.gtf \
--samples_file $proj_dir/Samples_Master_Tit-Reg-Pas_merge.txt \
--out_prefix DTU --aligner STAR > SuperTranscripts-DTU.log 2>&1


############################################
# 2.1. Evaluating the quality of the assembly
# Assembly metrics
############################################
#https://southgreenplatform.github.io/tutorials//bioanalysis/trinity/

perl $Trinity_path/util/align_and_estimate_abundance.pl \
--transcripts Trinity.fasta \
--seqType fq \
--samples_file $proj_dir/Samples_Master_Tit-Reg-Pas_merge.txt \
--est_method salmon \
--trinity_mode --prep_reference \
--output_dir outdir_estimate-ab > estimate-ab.log 2>&1

############################################
#DE_Trans

$Trinity_path/util/abundance_estimates_to_matrix.pl --est_method salmon \
--gene_trans_map Trinity.fasta.gene_trans_map \
--out_prefix Trinity_trans \
--name_sample_by_basedir \
Tit_rep1/quant.sf \
Tit_rep2/quant.sf \
Tit_rep3/quant.sf \
Tit_rep4/quant.sf \
Tit_rep5/quant.sf \
Tit_rep6/quant.sf \
Reg_rep1/quant.sf \
Reg_rep2/quant.sf \
Reg_rep3/quant.sf \
Reg_rep4/quant.sf \
Reg_rep5/quant.sf \
Reg_rep6/quant.sf \
Pas_rep1/quant.sf \
Pas_rep2/quant.sf \
Pas_rep3/quant.sf \
Pas_rep4/quant.sf \
Pas_rep5/quant.sf \
Pas_rep6/quant.sf > Trinity_trans_salmon.log 2>&1

############################################
#DE_gene


#Compute N50 based on the top-most highly expressed transcripts (Ex50)
$Trinity_path/util/misc/contig_ExN50_statistic.pl Trinity_trans.isoform.TMM.EXPR.matrix Trinity.fasta > ExN50.stats

#Tools to evaluate transcriptomes
$Trinity_path/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta > Trinity.longest.fasta

############################################
# 2.2. Identifying differentially expressed (DE) transcripts
############################################
#Extracting differentially expressed transcripts and generating heatmaps

############################################
$Trinity_path/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix Trinity_trans.isoform.counts.matrix \
--samples_file $proj_dir/Samples_Master_Tit-Reg-Pas_merge.txt \
--method DESeq2 \
--output DESeq2_trans > Trinity_DE-run.log 2>&1

#Extracting differentially expressed transcripts and generating heatmaps
#Extract those differentially expressed (DE) transcripts that are at least 4-fold (C is set to 2^(2) ) differentially expressed at a significance of <= 0.001 (-P 1e-3) in any of the pairwise sample comparisons

cd DESeq2_trans/
$Trinity_path/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../Trinity_trans.isoform.TMM.EXPR.matrix \
--samples $proj_dir/Samples_Master_Tit-Reg-Pas_merge.txt -P 1e-3 -C 2 > Trinity_DE-analyze.log 2>&1

############################################
#Run the DE analysis at the gene level

$Trinity_path/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix Trinity_trans.gene.counts.matrix \
--samples_file $proj_dir/Samples_Master_Tit-Reg-Pas_merge.txt \
--method DESeq2 \
--output DESeq2_genes > Trinity_DGE-run.log 2>&1

#Extracting differentially expressed transcripts and generating heatmaps
#Extract those differentially expressed (DE) transcripts that are at least 4-fold (C is set to 2^(2) ) differentially expressed at a significance of <= 0.001 (-P 1e-3) in any of the pairwise sample comparisons
cd DESeq2_genes/
$Trinity_path/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix ../Trinity_trans.gene.TMM.EXPR.matrix \
--samples $proj_dir/Samples_Master_Tit-Reg-Pas_merge.txt -P 1e-3 -C 2 > Trinity_DGE-analyze.log 2>&1

############################################


echo "Script done...."


