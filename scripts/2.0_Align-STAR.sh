#!/bin/bash


#Rna-seq using STAR
#./2.0_Align-STAR.sh > 2.0_Align-STAR.sh.log 2>&1


echo "run script for rna-seq-analysis"
##########################################
#STAR-Aligner
##########################################

#make REF index
<<COMMENT

################################################################################
########## Chenopodium_Quinoa

#http://malooflab.phytonetworks.org/wiki/bioinformatics/

ref=/home/gala0002/proj/proj_dir/REF_Genome/Ref_Chenopodium_Quinoa/

nice -n 5 STAR --runMode genomeGenerate \
--genomeDir ${ref} \
--genomeFastaFiles ${ref}GCF_001683475.1_ASM168347v1_genomic.fna \
--sjdbGTFfile ${ref}GCF_001683475.1_ASM168347v1_genomic.gff \
--runThreadN 60 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfeatureExon CDS \
--genomeSAindexNbases 12 \
--limitGenomeGenerateRAM 300000000000

##########

gffread GCF_001683475.1_ASM168347v1_genomic.gff -g GCF_001683475.1_ASM168347v1_genomic.fna -w GCF_001683475.1_ASM168347v1_trans.fa

##########

COMMENT

ulimit -n 65535
ulimit -c unlimited

#Make STAR alignment
######################
work_dir=/home/gala0002/proj/proj_dir/
ref=/home/gala0002/proj/proj_dir/REF_Genome/Ref_Chenopodium_Quinoa/

mkdir -p ${work_dir}NG-14833_2.0_Align-STAR_Quinoa_Genome/
out_dir=${work_dir}NG-14833_2.0_Align-STAR_Quinoa_Genome/


######################################################
#####---STAR with 2 Lanes:  P16061
######################################################

cd ${work_dir}NG-14833_1.1_sort-trim/

#for nbr in Sample_484-10-1/; do
#for nbr in Sample*/; do

for nbr in `ls /home/gala0002/proj/proj_dir/NG-14833_1.1_sort-trim/`

do
echo "Sample_DIR: $nbr"/

#nbr1=$(echo $nbr | sed 's=/[^/]*$==;s/\.$//')
nbr1=$(echo $nbr | cut -d"-" -f1)

echo $nbr1
#echo $nbr2

mkdir -p ${out_dir}${nbr1}/

temp_dir=${out_dir}${nbr1}/
cd ${temp_dir}

#All-reads
echo "Processing sample: ${nbr1}"

#5.2.2) Make START alignment

nice -n 5 STAR --genomeDir ${ref} \
--readFilesIn ${work_dir}NG-14833_1.1_sort-trim/${nbr}"/"${nbr1}-sort-trim_1.fq.gz ${work_dir}NG-14833_1.1_sort-trim/${nbr}"/"${nbr1}-sort-trim_2.fq.gz \
--sjdbGTFfile ${ref}GCF_001683475.1_ASM168347v1_genomic.gff \
--outFileNamePrefix ${temp_dir}${nbr1}-sort-trim-STAR \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts \
--twopassMode Basic \
--alignIntronMax 15000 \
--outFilterIntronMotifs RemoveNoncanonical \
--runThreadN 66 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfeatureExon CDS \
--outReadsUnmapped Fastx \
--readFilesCommand zcat

#nice -n 5 samtools sort -n ${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome.out.bam -o ${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome_sort.bam -@ 65

<<COMMENT
##nice -n 5 samtools sort -n ${temp_dir}${nbr1}-sort-trim-STARAligned.sortedByCoord.out.bam -o ${temp_dir}${nbr1}-sort-trim-STARAligned_sort.bam -@ 65
# remove all other files 

rm ${temp_dir}${nbr1}-sort-trim-STARAligned.sortedByCoord.out.bam
rm ${temp_dir}${nbr1}-sort-trim-STARLog.out
rm ${temp_dir}${nbr1}-sort-trim-STARLog.progress.out
rm ${temp_dir}${nbr1}-sort-trim-STARSJ.out.tab
rm -fr ${temp_dir}${nbr1}-sort-trim-STAR_STARgenome/
rm -fr ${temp_dir}${nbr1}-sort-trim-STAR_STARpass1/

# sorting mapped $ID reads vs $ID assembled contigs with 20 threads (throws error not working )
#samtools sort -n ${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome.out.bam -o ${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome_sort.bam -@ 65

# sort bam with picard 
java -jar /bioinfo/picard/picard.jar SortSam \
      I=${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome.out.bam \
      O=${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome_sort.bam \
      SORT_ORDER=coordinate

# index sorted bam with 20 threads
samtools index ${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome_sort.bam ${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome_sort.bam.bai -@ 65

# crop bam file 

samtools view -h ${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome_sort.bam rna31317 > /home/gala0002/proj/proj_dir/NG-14833_8.1_crop-bam/${nbr1}-sort-trim-STARAligned.toTranscriptome_sort_rna31317.bam
samtools view -h ${temp_dir}${nbr1}-sort-trim-STARAligned.toTranscriptome_sort.bam rna1790 > /home/gala0002/proj/proj_dir/NG-14833_8.1_crop-bam/${nbr1}-sort-trim-STARAligned.toTranscriptome_sort_rna1790.bam

COMMENT

done


echo "Done for independent samples..."


echo "Script done...."

