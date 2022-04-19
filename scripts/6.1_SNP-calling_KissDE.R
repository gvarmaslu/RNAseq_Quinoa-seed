#!/usr/bin/Rscript


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("kissDE")
#or
#install.packages("/Volumes/Seagate/softwares/kissDE_1.4.0.tar.gz",repos=NULL, type="source")


library(kissDE)

# dir1="/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_Tit/"
# dir2="/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_Reg/"
# dir3="/home/gala0002/proj/proj_dir/NG-14833_2.0_Kissplice_Pas/"

dir_inout="/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_2.0_Kissplice_all/"


# TEST

snp<-kissplice2counts(paste(dir_inout,"results_cond1_replicat1_R1_cond1_replicat1_R2_cond1_replicat2_R1_cond1_replicat2_R2_cond2_replicat1_R1_cond2_replicat1_R2_cond2_replicat2_R1_cond2_replicat2_R2_k41_coherents_type_0.fa",sep=""), pairedEnd=TRUE )
human_conditions<-c("c1","c1","c2","c2")
res<-diffExpressedVariants(snp, human_conditions, pvalue=1)
write.table(res$finalTable, file=paste(dir_inout,"kissDE_output",sep=""), sep="\t", quote=FALSE)


#NG-14833_Tit
snp<-kissplice2counts(paste(dir_inout,"results_NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_1_NG-14833_Tit_15dpa_1_lib233705_5747_1-sortmerna-trimmomatic_2_Merge_1_Merge_2_NG-14833_Tit_15dpa_3_lib233707_5747_1-sortmerna-trimmomatic_k41_coherents_type_0a.fa",sep=""), pairedEnd=TRUE )
human_conditions<-c("15dpa","15dpa","15dpa","25dpa","25dpa","25dpa")
res<-diffExpressedVariants(snp, human_conditions, pvalue=1)
write.table(res$finalTable, file=paste(dir_inout,"kissDE_Tit_output",sep=""), sep="\t", quote=FALSE)


#NG-14833_Reg
snp<-kissplice2counts(paste(dir_inout,"results_NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_1_NG-14833_Reg_15dpa_1_lib233717_5767_1-sortmerna-trimmomatic_2_Merge_1_Merge_2_Merge_1_Merge_2_Merge_1_Merge_2_NG-14833_Reg_25dpa_2_lib2337k41_coherents_type_0a.fa",sep=""), pairedEnd=TRUE )
human_conditions<-c("15dpa","15dpa","15dpa","25dpa","25dpa","25dpa")
res<-diffExpressedVariants(snp, human_conditions, pvalue=1)
write.table(res$finalTable, file=paste(dir_inout,"kissDE_Reg_output",sep=""), sep="\t", quote=FALSE)


#NG-14833_Pas
snp<-kissplice2counts(paste(dir_inout,"results_NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_1_NG-14833_Pas_E_1_lib233711_5747_1-sortmerna-trimmomatic_2_NG-14833_Pas_E_2_lib233712_5747_1-sortmerna-trimmomatic_1_NG-14833_Pas_E_2_lib233712k41_coherents_type_0a.fa",sep=""), pairedEnd=TRUE )
human_conditions<-c("15dpa","15dpa","15dpa","25dpa","25dpa","25dpa")
res<-diffExpressedVariants(snp, human_conditions, pvalue=1)
write.table(res$finalTable, file=paste(dir_inout,"kissDE_Pas_output",sep=""), sep="\t", quote=FALSE)


# script done...

