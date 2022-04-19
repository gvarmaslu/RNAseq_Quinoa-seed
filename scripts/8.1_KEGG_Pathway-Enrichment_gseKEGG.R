#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")
#install.packages("clusterProfiler")

######################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("KEGGprofile")
BiocManager::install("enrichKEGG")

library(clusterProfiler)
library(dplyr)
library(tidyverse)

library(enrichKEGG)
library(KEGGprofile)
library(KEGGREST)

###############

# gseMKEGG or enrichMKEGG

#https://yulab-smu.github.io/clusterProfiler-book/chapter6.html#kegg-gene-set-enrichment-analysis

#Chapter 6 KEGG analysis
#https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#######
#3.1 Input data
#######
#d <- read.csv(your_csv_file)
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/Reg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv")
#DE.table <- read.delim("/Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv")
DE.table <- read.delim("/Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars_GeneIDs_FC.txt")
#########################################################
#Prepare Input
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
# remove duplicates 
DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("AT_hit")]),]
DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$AT_hit == "No_AT_hit",]


DE.table_dedup_ids = DE.table[!duplicated(DE.table[c("Gene_ID")]),]
#DE.table_dedup_ids = DE.table_dedup_ids[!DE.table_dedup_ids$Gene_ID == ".",]

# Create a vector of the gene unuiverse
kegg_gene_list <- DE.table_dedup_ids$log2FoldChange

# Name vector with ENTREZ ids
#dedup_ids_ext_GID <- data.frame(do.call("rbind", strsplit(as.character(DE.table_dedup_ids$Gene_ID), ":", fixed = TRUE)))
names(kegg_gene_list) <- DE.table_dedup_ids$Gene_ID
names(kegg_gene_list) <- DE.table_dedup_ids$AT_hit

# omit any NA values 
kegg_gene_list <- na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

#http://yulab-smu.top/clusterProfiler-book/chapter6.html#kegg-over-representation-test
#6.1 KEGG over-representation test
# gene <-names(kegg_gene_list)[abs(kegg_gene_list) > 2]

#gene <- na.omit(dedup_ids_ext_GID$X2)

gene <- na.omit(DE.table_dedup_ids$AT_hit)

# quinoa: cqi 
# Arabidopsis thaliana: ath

# 6.1 KEGG over-representation test
kk <- enrichKEGG(gene         = gene,
                 organism     = 'ath',
                 pvalueCutoff = 0.05)
head(kk)


#Create gseKEGG object
#6.2 KEGG Gene Set Enrichment Analysis

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'ath',
               minGSSize    = 3,
               maxGSSize = 100,
               pvalueCutoff = 0.05, verbose      = FALSE)

head(kk2, n = 20L)

#####

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = 'ath',
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
head(kk2, n = 20L)

#####

#6.3 KEGG Module over-representation test

mkk <- enrichMKEGG(gene = gene,
                   organism = 'ath')

# 6.4 KEGG Module Gene Set Enrichment Analysis

mkk2 <- gseMKEGG(geneList = kegg_gene_list,
                 organism = 'ath')

#Merge DF

outdat_PM <- kk2
outdat_delta <- kk2

joined_df <- merge(outdat_PM, outdat_delta, by.x = "ID", 
                   by.y = "ID", all.x = T, all.y = T)

###plot

require(DOSE)
#dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
#dotplot(outdat_PM , showCategory = 24, title = "Enriched Pathways in PMTVWT" , split=".sign") + facet_grid(.~.sign)
#dotplot(outdat_delta, showCategory = 24, title = "Enriched Pathways in Delta8K" , split=".sign") + facet_grid(.~.sign)



#dataframe3 <- merge(outdat_PM, outdat_delta, by=c("pathway.code"),all=TRUE)

#write.table(kk, file = "/Volumes/Mac_HD2/proj_dir/RESULTS_2021/Reg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars_gseKEGG_join.tsv", row.names=T, sep="\t")

write.table(kk, file = "/Volumes/Mac_HD2/proj_dir/RESULTS_2021/Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars_gseKEGG_join.tsv", row.names=T, sep="\t")

#########################################################

DE.table <- read.delim("/Volumes/Mac_HD2/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv")

dedup_ids = DE.table[!duplicated(DE.table[c("Gene_ID")]),]

# Create a vector of the gene unuiverse
kegg_gene_list <- dedup_ids$log2FoldChange

# 
dedup_ids_ext_GID <- data.frame(do.call("rbind", strsplit(as.character(dedup_ids$Gene_ID), ":", fixed = TRUE)))

names(kegg_gene_list) <- dedup_ids_ext_GID$X2



## assume that 1st column is ID
## 2nd column is fold change
## feature 1: numeric vector
#geneList <- DE.table[,2]

## feature 2: named vector
#names(geneList) <- as.character(DE.table[,1])

## feature 3: decreasing order
geneList <- sort(kegg_gene_list, decreasing = TRUE)

#######

#6.1 KEGG over-representation test
#data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 1]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'cqi',
                 pvalueCutoff = 0.05)
head(kk)

### Merge

outdat_PM <- kk
outdat_delta <- kk

joined_df <- merge(outdat_PM, outdat_delta, by.x = "ID", 
                   by.y = "ID", all.x = T, all.y = T)

#dataframe3 <- merge(outdat_PM, outdat_delta, by=c("pathway.code"),all=TRUE)

write.table(joined_df, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.delta8K-PMTVWT-vs-Mock_enrichKEGG_join_all-0.05.tsv", row.names=T, sep="\t")


#6.3 KEGG Module over-representation test
mkk <- enrichMKEGG(gene = gene,
                   organism = 'cqi')

head(mkk)
#6.4 KEGG Module Gene Set Enrichment Analysis

mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'cqi')
head(mkk2)

#write.table(kk, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_PMTVWT_anno_pars_KEGG-In_enrichKEGG.tsv", row.names=T, sep="\t")
#write.table(kk, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.P5e-2_C1.0_Mock_vs_delta8K_anno_pars_KEGG-In_enrichKEGG.tsv", row.names=T, sep="\t")

write.table(kk, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.all_Mock_vs_PMTVWT_KEGG-In_enrichKEGG.tsv", row.names=T, sep="\t")
#write.table(kk, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.all_Mock_vs_delta8K_KEGG-In_enrichKEGG.tsv", row.names=T, sep="\t")


#########
#https://yulab-smu.github.io/clusterProfiler-book/chapter12.html#browsekegg
#Chapter 12 Visualization of Functional
#12.11 browseKEGG

browseKEGG(kk, 'ath03010')
browseKEGG(kk2, 'cqi00941')


browseKEGG(kk, 'ath00941')
#12.12 pathview from pathview package
library("pathview")

nta04712 <- pathview(gene.data  = geneList,
                     pathway.id = 'ath03010',
                     species    = "ath",
                     limit      = list(gene=max(abs(geneList)), cpd=1))


####################################################################################
#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#clusterprofiler

#KEGG enrichment analysis

#####
geneList <- DE.table %>% filter(P.Value < "0.05", !is.na(P.Value), abs(FC) < "1" ) 


## feature 1: numeric vector
genes <- geneList[,2]

## feature 2: named vector
names(genes) <- as.character(geneList[,1])

## feature 3: decreasing order
genes <- sort(genes, decreasing = TRUE)

##
kk <- enrichKEGG(gene = names(genes), organism = 'nta')
head(kk, n=10)



########
#https://www.infoworld.com/article/3454356/how-to-merge-data-in-r-using-r-merge-dplyr-or-datatable.html
#Merge DF

outdat_PM <- outdat
outdat_delta <- outdat

joined_df <- merge(outdat_PM, outdat_delta, by.x = "pathway.code", 
                   by.y = "pathway.code", all.x = T, all.y = T)

#dataframe3 <- merge(outdat_PM, outdat_delta, by=c("pathway.code"),all=TRUE)

write.table(joined_df, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v3/all_DESeq2.DE_results.all_Mock_vs_delta8K-PMTVWT_KEGG-In_KEGGREST_join_all.tsv", row.names=T, sep="\t")



#############################################

