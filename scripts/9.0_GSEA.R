#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")
#install.packages("clusterProfiler")

######################

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# #BiocManager::install("KEGGprofile")
# BiocManager::install("enrichKEGG")


##
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("org.At.tair.db")



library(clusterProfiler)
library(dplyr)
library(tidyverse)

#library(enrichKEGG)
library(KEGGprofile)
library(KEGGREST)
require(DOSE)
###############

# gseMKEGG or enrichMKEGG

#https://yulab-smu.github.io/clusterProfiler-book/chapter6.html#kegg-gene-set-enrichment-analysis

#Chapter 6 KEGG analysis
#https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#######
#3.1 Input data
#######
#d <- read.csv(your_csv_file)

#lst1=c("LfPo_vs_LfCl","LfWt_vs_LfCl","RtPo_vs_RtCl","RtWt_vs_RtCl")
#lst1=c("Po_vs_Cl","WtPo_vs_Cl","Wt_vs_Cl")
lst1=c("Pas_vs_Reg-Tit")
#lst1=c("Pas_vs_Reg")
#lst1=c("Pas_vs_Tit")

#for (i in lst1){
#  print(i)
i = lst1
filename=paste(lst1,".DESeq2.DE_results.P5e-2_C0.DE.subset_anno_AT-GIs.tsv" , sep="")
workdir="/Volumes/Mac_HD2/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/"
#workdir="/Volumes/Mac_HD2/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg/"
#workdir="/Volumes/Mac_HD2/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Tit/"
#workdir="/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_sugar-beet/work/DESeq2_genes_WtPo-Cl/"

DE.table <- read.delim(paste(workdir,filename, sep="/"))

# setting up the main directory
main_dir <- "/Volumes/Mac_HD2/proj_dir/KEGG_analysis_2022"
# setting up the sub directory
#sub_dir <- "LfPo_vs_LfCl_P5e-2_C0.DE"
sub_dir <- paste(i,"_P5e-2_C0.DE" , sep="")

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

#6.1 KEGG over-representation test
# 
kk <- enrichKEGG(gene         = gene,
                organism     = 'ath',
                pvalueCutoff = 0.05,
                keyType = "ncbi-geneid")
# 
# kk
# 
# browseKEGG(kk2, 'bvg01230')

#6.2 KEGG Gene Set Enrichment Analysis
#Create gseKEGG object
##https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#minGSSize    = 5,
#maxGSSize = 500,
#minGSSize    = 3,
#organism     = "cqi",
gse <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "ath",
               nPerm        = 10000,
               minGSSize    = 5,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",keyType = "ncbi-geneid")
#gse
#head(gse, n = 20L)
#browseKEGG(gse, 'bvg03040')

##### plot
#########################################################
# boxplot 
#dotplot(gse, showCategory = 20, title = "Enriched Pathways" , split=".sign", font.size = 16) + facet_grid(.~.sign) 

# save to file but check the format 
pdf(file="GSE_Pathways.pdf")
dotplot(gse, showCategory = 20, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
dev.off()

# Map
# Produce the native KEGG plot (PNG)
#dme <- pathview(gene.data=kegg_gene_list, pathway.id="ath01110", species = "bvg")
library(pathview)
#single sample 
#dme <- pathview(gene.data=kegg_gene_list, pathway.id="bvg01110", species = "bvg")

# PATHWAY VIEW 
# browseKEGG(gse, 'cqi00941')
#download.kegg(pathway.id="cqi00941", species = "cqi", kegg.dir = ".")

#dme <- pathview(gene.data=kegg_gene_list, pathway.id="cqi00941", species = "cqi", gene.idtype = "KEGG", kegg.native= T)
#gdata.osa=sim.mol.data(mol.type="gene", species ="cqi", id.type="kegg", nmol=10000)
#pv.out=pathview(gene.data = gdata.osa, pathway.id = "00195", gene.idtype = "KEGG", species = "Oryza sativa japonica")

# all samples 
for (i in gse@result$ID){
  ###
  print(i)
  dme <- pathview(gene.data=kegg_gene_list, pathway.id=i, species = "ath")
  ###
}

############
# remove xml files 
mydir <-file.path(main_dir, sub_dir)
delfiles <- dir(path=mydir,pattern="*xml")
file.remove(file.path(mydir, delfiles))

# Produce a different plot (PDF) (not displayed here)
#dme <- pathview(gene.data=kegg_gene_list, pathway.id="bvg03040", species = "bvg", kegg.native = F)
#knitr::include_graphics("bvg03040.pathview.png")

#########################################################
##### network analysis 
#library(DOSE)
#emapplot(gse, showCategory = 10)

### Category Netplot
library(ggnewscale)
pdf(file="GSE_Category-Netplot.pdf")
cnetplot(gse, categorySize="pvalue", foldChange=kegg_gene_list, showCategory = 6)
dev.off()

#cnetplot(gse, categorySize="pvalue", foldChange=kegg_gene_list, showCategory = 6)

#Ridgeplot
#Helpful to interpret up/down-regulated pathways.
library(ggridges)
#ridgeplot(gse) + labs(x="enrichment distribution")

#ridgeplot(gse)

pdf(file="up-down-regulated-pathways.pdf")
ridgeplot(gse)
dev.off()

# 

#########################################################

##### merging two results in a single file 
# outdat_PM <- kk2
# outdat_delta <- kk2
# joined_df <- merge(outdat_PM, outdat_delta, by.x = "ID", 
#                    by.y = "ID", all.x = T, all.y = T)
###plot
# require(DOSE)
# dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

# write result to file
#write.table(joined_df, file = "/Volumes/Mac_HD2/proj_Ramesh/RNA-seq_2020/RESULTS_2020_v4/all_DESeq2.DE_results.delta8K-PMTVWT-vs-Mock_gseKEGG_join_all-0.05.tsv", row.names=T, sep="\t")
write.table(gse, file = "./KEGG_GSEA.tsv", row.names=F, sep="\t")
write.table(kk, file = "./KEGG_ORA.tsv", row.names=F, sep="\t")

#########################################################

#}

#############################################

