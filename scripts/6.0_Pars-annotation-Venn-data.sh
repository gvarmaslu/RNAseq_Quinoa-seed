#!/bin/bash


#Pars DEG

echo "run script "

##########################################
##########################################
work_dir=/home/gala0002/proj/proj_dir/NG-14833_6.0_Venn-diagram/

cd $work_dir

###########
# extract mRNA IDs

#cp ../DESeq2_genes_Quinoa_*/*.P1e-2_C1.DE.subset_anno.tsv ./

cat Pas_vs_Reg.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | cut -f1 | grep -v "Quinoa_mRNA-ID" > Pas_vs_Reg_P1e-2_C1_mRNA-IDs.tsv
cat Pas_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | cut -f1 | grep -v "Quinoa_mRNA-ID" > Pas_vs_Tit_P1e-2_C1_mRNA-IDs.tsv
cat Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | cut -f1 | grep -v "Quinoa_mRNA-ID" > Reg_vs_Tit_P1e-2_C1_mRNA-IDs.tsv

cat Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | cut -f1 | grep -v "Quinoa_mRNA-ID" > Pas_vs_Reg-Tit_P1e-2_C1_mRNA-IDs.tsv
cat AReg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | cut -f1 | grep -v "Quinoa_mRNA-ID" > Reg_vs_Pas-Tit_P1e-2_C1_mRNA-IDs.tsv
cat ATit_vs_Pas-Reg.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | cut -f1 | grep -v "Quinoa_mRNA-ID" > Tit_vs_Pas-Reg_P1e-2_C1_mRNA-IDs.tsv

###########
# Extract all gene IDs

cat Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | cut -f33 | grep -w -v "Gene_Name" | grep -w -v "." | sort -u > Pas_vs_Reg-Tit_P1e-2_C1_Gene-IDs.tsv
cat AReg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | cut -f33 | grep -w -v "Gene_Name" | grep -w -v "." | sort -u > Reg_vs_Pas-Tit_P1e-2_C1_Gene-IDs.tsv
cat ATit_vs_Pas-Reg.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | cut -f33 | grep -w -v "Gene_Name" | grep -w -v "." | sort -u > Tit_vs_Pas-Reg_P1e-2_C1_Gene-IDs.tsv

# Extract all unique gene IDs, remove empty col (.) and with full content

#cat Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv | awk -F"\t" '!_[$33]++' | awk '$31 != "."' > venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvReg_1371_anno_uniq.tsv

###########
### Filter annotation files 
work_dir=/Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/

cd $work_dir

#cp ../NG-14833_6.0_Venn-diagram/*.P1e-2_C1.DE.subset_anno.tsv ./
# rename 
cat AReg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | sed 's/AReg//g' > Reg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv
cat ATit_vs_Pas-Reg.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | sed 's/ATit//g' > Tit_vs_Pas-Reg.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv

# remove 
rm AReg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv
rm ATit_vs_Pas-Reg.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv




#Pas_vs_Reg-Tit Reg_vs_Pas-Ti Tit_vs_Pas-Reg Pas_vs_Reg Pas_vs_Tit Reg_vs_Tit
#Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg 
#Pas_vs_Reg Pas_vs_Tit Reg_vs_Tit Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg 

for i in Pas_vs_Reg Pas_vs_Tit Reg_vs_Tit Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg

do 
echo $i 

# take unique genes and remove missing genes 

#cat Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | awk -F"\t" '!_[$33]++' | awk '$33 != "."' > Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq.tsv
cat Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-10_C2.DE.subset_anno.tsv | awk -F"\t" '!_[$33]++' | awk '$33 != "."' > Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-10_C2.DE.subset_anno_uniq.tsv

#cat ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv | awk -F"\t" '!_[$33]++' | awk '$33 != "."' > ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq.tsv 

################################################
#Pars with Python script

#python 7.3_Pars-annotation-Venn-data-pars.py Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq.tsv Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv /Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/
python /Volumes/Mac_HD2/RNAseq-analysis-run_Backup/scripts_quinoa_v2/7.3_Pars-annotation-Venn-data-pars.py ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq.tsv ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv /Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/  

# Pars : count TF database genes
echo "total:"
cat ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv | wc -l 
echo "PlantTFDBv5.0:"
cat ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv | cut -f49-51 | grep "|" | wc -l 
echo "ITAKv18.12:"
cat ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv | cut -f52-55 | grep "|" | wc -l 
echo "InAT:"
cat ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv | cut -f44 | awk '$1 !="."' | wc -l 
cat ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv | cut -f59 | awk '$1 !="No_AT_hit"' | wc -l 

done	

echo "Script done...."


#################################################
##########################
#SNP calling filter on Genes (whole genome SNP calling)
for i in Pas_vs_Reg Pas_vs_Tit Reg_vs_Tit Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg
do 
echo $i 
cat ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq.tsv | cut -f33 | grep -v "Gene_Name" > ${i}_Genes_uniq.tsv
done

#SNP calling filter on mRNA (whole genome SNP calling)
for i in Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg
do 
echo $i 
cat ${i}.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq.tsv | cut -f1 | grep -v "Quinoa_mRNA-ID" > ${i}_mRNA_uniq.tsv
done

# total SNP count 
#for i in Pas Reg Tit; do  echo $i; cat ${SNP_DIR}variants_${i}/filtered_output_sID_flr2_snpeff_Prot.vcf | grep -v "^#" | wc -l; done

#Tit_vs_Pas-Reg

for i in Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg
do 
echo $i
### split string

jj=$(echo $i | cut -d"_" -f1)
echo $jj

##### 
#SNP_DIR=/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_2.1_STAR_Genome-CDS-ass/
SNP_DIR=/mnt/udda-backup/data/RAW-DATA-BACKUP/Asa_bu/NG-14833_2.1_STAR_Genome-ass/

grep -w -f ${i}_Genes_uniq.tsv ${SNP_DIR}variants_${jj}/filtered_output_sID_flr2_snpeff_Prot.vcf | grep -v "#" > ${i}_Genes_uniq_SNPsel_${jj}.tsv

###### make VCF 
cat ${SNP_DIR}variants_${jj}/filtered_output_sID_flr2_snpeff_Prot.vcf | grep "^#" > head_${jj}

## merge vcf header
cat head_${jj} ${i}_Genes_uniq_SNPsel_${jj}.tsv > ${i}_Genes_uniq_SNPsel_${jj}.vcf


##
rm head_*

done

#######################################################################################
####### mRNA selection

for i in Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg
do 
echo $i
### split string

jj=$(echo $i | cut -d"_" -f1)
echo $jj

##### 
#SNP_DIR=/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_2.1_STAR_Genome-CDS-ass/
SNP_DIR=/mnt/udda-backup/data/RAW-DATA-BACKUP/Asa_bu/NG-14833_2.1_STAR_Genome-ass/

<<COMMENT
####### Filter VCF file using SNPEFF
#Filter alignment artifacts to get promising SNP
#SNPSIFT filter  
cat ${SNP_DIR}variants_${jj}/filtered_output_sID_flr2_snpeff_Prot.vcf | java -jar /bioinfo/SNPEFF/snpEff/SnpSift.jar filter "(( QUAL >= 30 ) && ( DP >= 10 )) && ( GEN[*].GQ >= 20)" > ${jj}_filtered_output_sID_flr2_snpeff_Prot_flr3.vcf

### Pars futher on VCF file 
cat ${jj}_filtered_output_sID_flr2_snpeff_Prot_flr3.vcf | grep -v "#" | awk -v FS='\t' -v OFS='\t' '{ split($8,a,"|"); print $1,$2,$3,$4,$5,$6,$7,$10,a[2],a[7]}' > ${i}_mRNA_uniq_SNPsel_${jj}.tsv

### Find list of genes in SNP file 

grep -w -f ${i}_mRNA_uniq.tsv ${i}_mRNA_uniq_SNPsel_${jj}.tsv > ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE.tsv

### find High-impact variants 
grep -w -f HIGH-MODERATE_impact_variants.txt ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE.tsv > ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv 

#cat ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_High-impact.tsv | awk -v FS="\t" -v OFS="\t" '{ split($8,a,":"); if((a[3]>=10) && ($6 >=30)) print $0; }' > ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_High-impact_qualflt.tsv
COMMENT

### extract filtered SNPs and make VCF file 

#grep -e "NW_018742244.1.*762784" Reg_filtered_output_sID_flr2_snpeff_Prot_flr3.vcf
cat ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv  | awk -v FS='\t' '{print $1".*"$2}' > ${jj}_search-list-pos.txt

grep -w -E -f ${jj}_search-list-pos.txt ${jj}_filtered_output_sID_flr2_snpeff_Prot_flr3.vcf > ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact_full.tsv

###### make VCF 
cat ${jj}_filtered_output_sID_flr2_snpeff_Prot_flr3.vcf | grep "^#" > head_${jj}

## merge vcf header
cat head_${jj} ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact_full.tsv > ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact_full.vcf


##
rm head_*

done

######

############################
######
#for i in Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg

for i in Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg
do
echo $i
#jj= $i | cut -d"_" -f1 
jj=$(echo $i | cut -d"_" -f1)
echo $jj
#echo "$i splice_acceptor_variant: "
grep -w "splice_acceptor_variant" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l; 
#echo "splice_donor_variant:"  
grep -w "splice_donor_variant" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l;
#echo "stop_gained:" 
grep -w "stop_gained" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l;
#echo "frameshift_variant:"
grep -w "frameshift_variant" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l;
#echo "stop_lost:" 
grep -w "stop_lost" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l; 
#echo "start_lost:"; 
grep -w "start_lost" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l; 

#echo "transcript_ablation:"; 
#grep -w "transcript_ablation" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l; 
#echo "inframe_insertion:"; 
#grep -w "inframe_insertion" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l; 
#echo "inframe_deletion:"; 
#grep -w "inframe_deletion" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l; 
#echo "missense_variant:"; 
grep -w "missense_variant" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l; 
#echo "protein_altering_variant:"; 
#grep -w "protein_altering_variant" ${i}_mRNA_uniq_SNPsel_${jj}_P1e-2_C1.DE_HIGH-MODERATE-impact.tsv | wc -l; 


#grep -f High_impact_variants.txt -w ${i}_Genes_uniq_SNPsel_${jj}.vcf > ${i}_Genes_uniq_SNPsel_High-impact_${jj}.tsv 
#cat ${i}_Genes_uniq_SNPsel_${jj}.vcf | grep "^#" > head_${jj}

## merge vcf header
#cat head_${jj} ${i}_Genes_uniq_SNPsel_High-impact_${jj}.tsv > ${i}_Genes_uniq_SNPsel_High-impact_${jj}.vcf

##
#rm head_*

done

###### 
for i in Pas_vs_Reg-Tit Reg_vs_Pas-Tit Tit_vs_Pas-Reg
do 
echo $i
#jj= $i | cut -d"_" -f1 
jj=$(echo $i | cut -d"_" -f1)
echo $jj

python 7.4_Pars-SNP-genes.py ${i}_Genes_uniq_SNPsel_High-impact_${jj}.tsv ${i}_Genes_uniq_SNPsel_High-impact_${jj}_pars.tsv /home/gala0002/proj/proj_dir/NG-14833_7.0_filter/

done
######

echo "done..."

#| awk '$6 >= 100'
###### 
for i in Pas_vs_Reg-Tit; 
do  echo $i; 
jj=$(echo $i | cut -d"_" -f1); 
echo $jj;
grep -w "frameshift_variant" ${i}_Genes_uniq_SNPsel_${jj}.vcf ;  

#grep -w "frameshift_variant" ${i}_Genes_uniq_SNPsel_${jj}.vcf | awk -F "\t" '{ print $1,$2,$3,$4,$5,$6,$7,$10,$2-100,$2+100 }';  
#grep -w "frameshift_variant" ${i}_Genes_uniq_SNPsel_${jj}.vcf | awk -F "\t" '{ split($8,a,"|"); print $1,$2,$3,$4,$5,$6,$7,$10,a[2],a[7]}';  

done


###### Crop bam files 

######

#samtools view -h example.bam 17:7512445-7513455
#NW_018745684.1:2621469-2624421

for i in NW_018745684.1:2621469-2624421;
do 
samtools view -h /mnt/udda-backup/data/RAW-DATA-BACKUP/Asa_bu/NG-14833_2.1_STAR_Genome-ass/variants_Pas/splitNCigar.bam ${i} > Pas_splitNCigar_${i}.sam
samtools view -h /mnt/udda-backup/data/RAW-DATA-BACKUP/Asa_bu/NG-14833_2.1_STAR_Genome-ass/variants_Reg/splitNCigar.bam ${i} > Reg_splitNCigar_${i}.sam
samtools view -h /mnt/udda-backup/data/RAW-DATA-BACKUP/Asa_bu/NG-14833_2.1_STAR_Genome-ass/variants_Tit/splitNCigar.bam ${i} > Tit_splitNCigar_${i}.sam

done 

###### extract DET geneIds from High-impact SNP ids 

#cat Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv | cut -f1 > Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars_rnaIDs.tsv
#grep -f Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars_rnaIDs.tsv -w Pas_vs_Reg-Tit_Genes_uniq_SNPsel_High-impact_Pas_pars.tsv > Pas_vs_Reg-Tit_Genes_uniq_SNPsel_High-impact_Pas_pars_P1e-2_C1.DE.tsv

#####
#for i in Pas_vs_Reg-Tit;  do  echo $i;  jj=$(echo $i | cut -d"_" -f1);  echo $jj;  grep -w "frameshift_variant" ${i}_Genes_uniq_SNPsel_${jj}.tsv | awk -F "\t" '{ split($8,a,"|"); print $1,$2,$3,$4,$5,$6,$7,$10,a[2],a[7]}';   done > Pas_vs_Reg-Tit_Genes_uniq_SNPsel_frameshift_variant.tsv

for i in Pas_vs_Reg-Tit;  do  echo $i;  jj=$(echo $i | cut -d"_" -f1);  echo $jj; cat ${i}_Genes_uniq_SNPsel_High-impact_${jj}.tsv | awk -v FS='\t' -v OFS='\t' '{ split($8,a,"|"); print $1,$2,$3,$4,$5,$6,$7,$10,a[2],a[7]}'; done > Pas_vs_Reg-Tit_Genes_uniq_SNPsel_High-impact_Pas_parsall.tsv

grep -f Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars_rnaIDs.tsv -w Pas_vs_Reg-Tit_Genes_uniq_SNPsel_High-impact_Pas_parsall.tsv > Pas_vs_Reg-Tit_Genes_uniq_SNPsel_High-impact_Pas_parsall_P1e-2_C1.DE.tsv


##### select DP of reads greater than 10 

#head Pas_vs_Reg-Tit_Genes_uniq_SNPsel_High-impact_Pas_parsall_P1e-2_C1.DE.tsv | awk -v FS="\t" -v OFS="\t" '{ split($8,a,":"); if(a[3]>=10); print $0 }' 

cat Pas_vs_Reg-Tit_Genes_uniq_SNPsel_High-impact_Pas_parsall_P1e-2_C1.DE.tsv | awk -v FS="\t" -v OFS="\t" '{ split($8,a,":"); if((a[3]>=10) && ($6 >=30)) print $0; }' > Pas_vs_Reg-Tit_Genes_uniq_SNPsel_High-impact_Pas_parsall_P1e-2_C1.DE_qualflt.tsv

grep -f interest_selIDs -w Pas_vs_Reg-Tit_Genes_uniq_SNPsel_High-impact_Pas_parsall_P1e-2_C1.DE_qualflt.tsv


######################################################
