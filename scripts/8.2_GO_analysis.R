#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")
#install.packages("clusterProfiler")

######################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("KEGGprofile")
BiocManager::install("enrichKEGG")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.At.tair.db")

library(DOSE)
library(AnnotationDbi)
library(org.At.tair.db)

library(clusterProfiler)
library(dplyr)
library(tidyverse)

library(enrichKEGG)
library(KEGGprofile)
library(KEGGREST)

###############
#https://shiring.github.io/rna-seq/deseq2/teaching/2016/09/29/DESeq2-course

#######
#3.1 Input data
#######
DE.table <- read.delim("/Volumes/Mac_HD2/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C0.5.DE.subset_anno_AT.tsv")

# remove duplicates 
DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("AT_hit")]),]

DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$AT_hit == "No_AT_hit",]

# split GeneIDs
#DE.table_dedup_ids$Gene_ID <- data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$Gene_ID), ":", fixed = TRUE)))$X2


#########################################################
#Prepare Input

geneList <- as.vector(DE.table_dedup_ids$log2FoldChange)
names(geneList) <- DE.table_dedup_ids$AT_hit
gene <- na.omit(DE.table_dedup_ids$AT_hit)


library(clusterProfiler)
#data(geneList, package="DOSE")
#gene <- names(geneList)[abs(geneList) > 2]
#gene.df <- bitr(gene, fromType = "TAIR",
#                toType = c("ENTREZID", "GO"),
#                OrgDb = org.At.tair.db)
#head(gene.df)

# Group GO

ggo <- clusterProfiler::groupGO(gene     = gene,
                                OrgDb    = org.At.tair.db,
                                ont      = "BP",
                                level    = 3)
head(summary(ggo)[,-5])

# GO over-representation test
ego <- enrichGO(gene          = geneList,
                OrgDb         = org.At.tair.db,
                keyType= 'TAIR',
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego)[,-8])


  ego <- enrichGO(gene = geneList,
                OrgDb = org.At.tair.db,
                ont = "CC",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.1,
                minGSSize = 10,
                maxGSSize = 500)




