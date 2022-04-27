#!/usr/bin/RScript


#https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html


library(clusterProfiler)
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

# Entrez gene ID
head(gene)

####
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("org.At.tair.db")

####
ggo <- groupGO(gene     = gene,
               OrgDb    = org.At.tair.db,
               keyType= 'TAIR',
               ont      = "CC")

head(ggo)

#https://github.com/YuLab-SMU/clusterProfiler/issues/282

library(org.At.tair.db)

load('/home/gala0002/proj/proj_dir/KEGG-analysis_2022/kall.rdata')

kallGOBP <- compareCluster(geneCluster = kall[1:2],
                           fun = 'enrichGO',
                           OrgDb = 'org.At.tair.db',
                           keyType= 'TAIR',
                           ont = 'BP',
                           universe = keys(org.At.tair.db),
                           pAdjustMethod = 'BH',
                           pvalueCutoff=0.05,
                           qvalueCutoff=0.1)

########

kallGOBP_2 <- enrichGO(kall[[1]],
         OrgDb = 'org.At.tair.db',
         keyType= 'TAIR',
         ont = 'BP',
         universe = keys(org.At.tair.db),
         pAdjustMethod = 'BH',
         pvalueCutoff=0.05,
         qvalueCutoff=0.1)

########
########



library(clusterProfiler)
library(dplyr)
library(tidyverse)

#library(enrichKEGG)
library(KEGGprofile)
library(KEGGREST)
require(DOSE)


#######
lst1=c("Pas_vs_Reg-Tit")
#lst1=c("Pas_vs_Reg")
#lst1=c("Pas_vs_Tit")

#######
i = lst1
filename=paste(lst1,".DESeq2.DE_results.P1e-2_C1.DE.subset_anno_AT.tsv" , sep="")
workdir="/home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/"
#workdir="/home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg/"
#workdir="/home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Tit/"

DE.table <- read.delim(paste(workdir,filename, sep="/"))

# setting up the main directory
main_dir <- "/home/gala0002/proj/proj_dir/GO_analysis_2022"
# setting up the sub directory
#sub_dir <- "LfPo_vs_LfCl_P5e-2_C0.DE"
sub_dir <- paste(i,"_P1e-2_C1.DE" , sep="")

# check if sub directory exists 
if (file.exists(sub_dir)){
  # specifying the working directory
  setwd(file.path(main_dir, sub_dir))
} else {
  # create a new sub directory inside
  # the main path
  dir.create(file.path(main_dir, sub_dir))
  # specifying the working directory
  setwd(file.path(main_dir, sub_dir))
}

#########################################################
#Prepare Input
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = DE.table[!duplicated(DE.table[c("Gene_ID")]),]

# Create a vector of the gene unuiverse
#kegg_gene_list <- dedup_ids$log2FoldChange

# Name vector with ENTREZ ids
#names(kegg_gene_list) <- dedup_ids$GeneID

###############
## cbind 
## match with the other list file

# wheat_data = data_frame(dedup_ids)
# wheat_bvg = data_frame(DE.table_bvg)
# df_flt <- merge(wheat_data, wheat_bvg, by = "GeneID_iwgsc_refseqv2.1")
# 
# glimpse(wheat_data)
# glimpse(wheat_bvg)
# glimpse(df_flt)

###############
# load complete file
df_flt = data_frame(dedup_ids)

# load to kegg variable
kegg_gene_list <- df_flt$log2FoldChange

# add KEGG IDs of CQI= Quinoa
#names(kegg_gene_list) <- df_flt$GeneID

#dedup_ids_ext_GID <- data.frame(do.call("rbind", strsplit(as.character(dedup_ids$Gene_ID), ":", fixed = TRUE)))
#names(kegg_gene_list) <- dedup_ids_ext_GID$X2

# add KEGG IDs of AT
#dedup_ids$Target_AT_PID_NCBIBLAST
#dedup_ids_ext_GID <-  df_flt$Target_AT_PID_NCBIBLAST
#dedup_ids_ext_GID_2 <- as.numeric(ifelse(dedup_ids_ext_GID==".", NA, dedup_ids_ext_GID))

df_flt$AT_hit[df_flt$AT_hit == "No_AT_hit"] <- NA
names(kegg_gene_list) <- df_flt$AT_hit

###############
# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list <- kegg_gene_list[!is.na(names(kegg_gene_list))]

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
gene <- names(kegg_gene_list)

###############
#GO-ORA
# pAdjustMethod = 'BH',
#minGSSize = 200,
kallGOBP <- enrichGO(gene,
                       OrgDb = 'org.At.tair.db',
                       keyType= 'TAIR',
                       ont = 'ALL',
                       universe = keys(org.At.tair.db),
                       pvalueCutoff=0.01,
                       qvalueCutoff=0.01)
dim(kallGOBP)
#6.4 GO Gene Set Enrichment Analysis
# GO-GSEA
#minGSSize    = 100,
#maxGSSize    = 500,
ego3 <- gseGO(geneList     = kegg_gene_list,
              OrgDb        = org.At.tair.db,
              keyType= 'TAIR',
              nPermSimple = 10000,
              ont          = "ALL",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

dim(ego3)
#
#6.6 Visualize enriched GO terms as a directed acyclic graph
#goplot(ego3)

#barplot(kallGOBP)
#clusterProfiler::dotplot(ego3, showCategory=10)
###############

#########################################################
# write result to file
write.table(kallGOBP, file = "./GO_ORA.tsv", row.names=F, sep="\t")
write.table(ego3, file = "./GO_GSEA.tsv", row.names=F, sep="\t")

###############
