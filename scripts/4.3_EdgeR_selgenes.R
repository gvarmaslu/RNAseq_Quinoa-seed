##########################
#EdgeR script
#Source protocol
#http://www.nathalievilla.org/doc/html/solution_edgeR-tomato-withcode.html

library("edgeR")
library(limma)
library(RColorBrewer)
library(mixOmics)

#source("http://bioconductor.org/biocLite.R")
#biocLite("HTSFilter", type = "source")

#BiocManager::install(version = "3.9")
#BiocManager::install("HTSFilter")
#BiocManager::install("mixOmics")

#4.0. Starting from count table
#Preparing the files (count table and design)

directory <- getwd()
dir(directory)

#Importing the files
# outdir_plot="/Volumes/Seagate/Backup-MAC_HD2/proj_Aakash/QE-1467_4.1_EdgeR_results_plot"
# outdir_work="/Volumes/Seagate/Backup-MAC_HD2/proj_Aakash/QE-1467_4.0_EdgeR_results/"

outdir_plot="/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_4.3_EdgeR_results_plot"
outdir_work="/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_4.3_EdgeR_results/"


dir.create(file.path(outdir_plot), showWarnings = F, recursive = FALSE)
dir.create(file.path(outdir_work), showWarnings = F, recursive = FALSE)
# 
# rawCountTable <- read.table("/Volumes/Seagate/Backup-MAC_HD2/proj_Aakash/QE-1467_3.2_express_combine-all/Express_counts_all_484-485-506-565-Nimbus-Stigg.csv", header=TRUE, sep=",", row.names=1, check.names=FALSE)
# sampleInfo <- read.table("/Volumes/Seagate/Backup-MAC_HD2/proj_Aakash/QE-1467_3.2_express_combine-all/Express_counts_all_484-485-506-565-Nimbus-Stigg_metadata_all.txt", header=TRUE, sep="\t", row.names=1)

rawCountTable <- read.table("/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_4.2_express_combine-all/Express_counts_all_Pas-Reg-Tit.csv", header=TRUE, sep=",", row.names=1, check.names=FALSE)
sampleInfo <- read.table("/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_4.2_express_combine-all/Express_counts_all_Pas-Reg-Tit_metadata_all.txt", header=TRUE, sep="\t", row.names=1)

#head(rawCountTable)
#head(sampleInfo)

################################ 
# #### Combinations 10d_R-vs-S, 12d_R-vs-S, R_10d-vs-12d, S_10d-vs-12d
# sampleInfo_10d_RvS <- sampleInfo[ which(sampleInfo$days=='10' ), ]
# sampleInfo_12d_RvS <- sampleInfo[ which(sampleInfo$days=='12'), ]
# sampleInfo_R_10dv12d <- sampleInfo[ which(sampleInfo$resp=='R'), ]
# sampleInfo_S_10dv12d <- sampleInfo[ which(sampleInfo$resp=='S'), ]
# 
# ###---select only specific samples; 10d_R-vs-S, 12d_R-vs-S, R_10d-vs-12d, S_10d-vs-12d
# rawCountTable_10d_RvS <- rawCountTable[ rownames(sampleInfo_10d_RvS) ]
# rawCountTable_12d_RvS <- rawCountTable[ rownames(sampleInfo_12d_RvS) ]
# rawCountTable_R_10dv12d <- rawCountTable[ rownames(sampleInfo_R_10dv12d) ]
# rawCountTable_S_10dv12d <- rawCountTable[ rownames(sampleInfo_S_10dv12d) ]
# 
# #4.1. Starting from count table
# dgeFull_10d_RvS <- DGEList(rawCountTable_10d_RvS, group=sampleInfo_10d_RvS$resp)
# dgeFull_12d_RvS <- DGEList(rawCountTable_12d_RvS, group=sampleInfo_12d_RvS$resp)
# dgeFull_R_10dv12d <- DGEList(rawCountTable_R_10dv12d, group=sampleInfo_R_10dv12d$days)
# dgeFull_S_10dv12d <- DGEList(rawCountTable_S_10dv12d, group=sampleInfo_S_10dv12d$days)
# 
# #4.1. Starting from count table ; Combinations 10d-vs-12d, R-vs-S, P-vs-O
# dgeFull_10dv12d <- DGEList(rawCountTable, group=sampleInfo$days)
# dgeFull_RvS <- DGEList(rawCountTable, group=sampleInfo$resp)
# dgeFull_PvO <- DGEList(rawCountTable, group=sampleInfo$po)

################################  

# #### Combinations "Pas_15dv25d","Reg_15dv25d","Tit_15dv25d"
# #### Combinations "Pas_15dvReg_15d","Pas_15dvTit_15d","Reg_15dvTit_15d","Pas_25dvReg_25d","Pas_25dvTit_25d","Reg_25dvTit_25d"
# #### Combinations "PasvReg","PasvTit","RegvTit"

sampleInfo_Pas_15dv25d <- sampleInfo[ which(sampleInfo$sampID=='Pas'), ]
sampleInfo_Reg_15dv25d <- sampleInfo[ which(sampleInfo$sampID=='Reg'), ]
sampleInfo_Tit_15dv25d <- sampleInfo[ which(sampleInfo$sampID=='Tit'), ]

sampleInfo_Pas_15dvReg_15d <- subset(sampleInfo, (days=='15dpa') & (sampID=='Pas' | sampID=='Reg'))
sampleInfo_Pas_15dvTit_15d <- subset(sampleInfo, (days=='15dpa') & (sampID=='Pas' | sampID=='Tit'))
sampleInfo_Reg_15dvTit_15d <- subset(sampleInfo, (days=='15dpa') & (sampID=='Reg' | sampID=='Tit'))

sampleInfo_Pas_25dvReg_25d <- subset(sampleInfo, (days=='25dpa') & (sampID=='Pas' | sampID=='Reg'))
sampleInfo_Pas_25dvTit_25d <- subset(sampleInfo, (days=='25dpa') & (sampID=='Pas' | sampID=='Tit'))
sampleInfo_Reg_25dvTit_25d <- subset(sampleInfo, (days=='25dpa') & (sampID=='Reg' | sampID=='Tit'))

sampleInfo_PasvReg <- sampleInfo[ which(sampleInfo$sampID=='Pas' | sampleInfo$sampID=='Reg'), ]
sampleInfo_PasvTit <- sampleInfo[ which(sampleInfo$sampID=='Pas' | sampleInfo$sampID=='Tit'), ]
sampleInfo_RegvTit <- sampleInfo[ which(sampleInfo$sampID=='Reg' | sampleInfo$sampID=='Tit'), ]

# ###---select only specific samples; "Pas_15dv25d","Reg_15dv25d","Tit_15dv25d"
# ###---select only specific samples; "Pas_15dvReg_15d","Pas_15dvTit_15d","Reg_15dvTit_15d","Pas_25dvReg_25d","Pas_25dvTit_25d","Reg_25dvTit_25d"
# ###---select only specific samples; "PasvReg","PasvTit","RegvTit"

rawCountTable_Pas_15dv25d <- rawCountTable[ rownames(sampleInfo_Pas_15dv25d) ]
rawCountTable_Reg_15dv25d <- rawCountTable[ rownames(sampleInfo_Reg_15dv25d) ]
rawCountTable_Tit_15dv25d <- rawCountTable[ rownames(sampleInfo_Tit_15dv25d) ]

rawCountTable_Pas_15dvReg_15d <- rawCountTable[ rownames(sampleInfo_Pas_15dvReg_15d) ]
rawCountTable_Pas_15dvTit_15d <- rawCountTable[ rownames(sampleInfo_Pas_15dvTit_15d) ]
rawCountTable_Reg_15dvTit_15d <- rawCountTable[ rownames(sampleInfo_Reg_15dvTit_15d) ]

rawCountTable_Pas_25dvReg_25d <- rawCountTable[ rownames(sampleInfo_Pas_25dvReg_25d) ]
rawCountTable_Pas_25dvTit_25d <- rawCountTable[ rownames(sampleInfo_Pas_25dvTit_25d) ]
rawCountTable_Reg_25dvTit_25d <- rawCountTable[ rownames(sampleInfo_Reg_25dvTit_25d) ]

rawCountTable_PasvReg <- rawCountTable[ rownames(sampleInfo_PasvReg) ]
rawCountTable_PasvTit <- rawCountTable[ rownames(sampleInfo_PasvTit) ]
rawCountTable_RegvTit <- rawCountTable[ rownames(sampleInfo_RegvTit) ]

#4.1. Starting from count table;  "Pas_15dv25d","Reg_15dv25d","Tit_15dv25d"
#4.1. Starting from count table;  "Pas_15dvReg_15d","Pas_15dvTit_15d","Reg_15dvTit_15d","Pas_25dvReg_25d","Pas_25dvTit_25d","Reg_25dvTit_25d"
dgeFull_Pas_15dv25d <- DGEList(rawCountTable_Pas_15dv25d, group=sampleInfo_Pas_15dv25d$days)
dgeFull_Reg_15dv25d <- DGEList(rawCountTable_Reg_15dv25d, group=sampleInfo_Reg_15dv25d$days)
dgeFull_Tit_15dv25d <- DGEList(rawCountTable_Tit_15dv25d, group=sampleInfo_Tit_15dv25d$days)

dgeFull_Pas_15dvReg_15d <- DGEList(rawCountTable_Pas_15dvReg_15d, group=sampleInfo_Pas_15dvReg_15d$sampID)
dgeFull_Pas_15dvTit_15d <- DGEList(rawCountTable_Pas_15dvTit_15d, group=sampleInfo_Pas_15dvTit_15d$sampID)
dgeFull_Reg_15dvTit_15d <- DGEList(rawCountTable_Reg_15dvTit_15d, group=sampleInfo_Reg_15dvTit_15d$sampID)

dgeFull_Pas_25dvReg_25d <- DGEList(rawCountTable_Pas_25dvReg_25d, group=sampleInfo_Pas_25dvReg_25d$sampID)
dgeFull_Pas_25dvTit_25d <- DGEList(rawCountTable_Pas_25dvTit_25d, group=sampleInfo_Pas_25dvTit_25d$sampID)
dgeFull_Reg_25dvTit_25d <- DGEList(rawCountTable_Reg_25dvTit_25d, group=sampleInfo_Reg_25dvTit_25d$sampID)


#4.1. Starting from count table ; Combinations "15dv25d"
#4.1. Starting from count table ; #"PasvReg","PasvTit","RegvTit"
dgeFull_15dv25d <- DGEList(rawCountTable, group=sampleInfo$days)

dgeFull_PasvReg <- DGEList(rawCountTable_PasvReg, group=sampleInfo_PasvReg$sampID )
dgeFull_PasvTit <- DGEList(rawCountTable_PasvTit, group=sampleInfo_PasvTit$sampID)
dgeFull_RegvTit <- DGEList(rawCountTable_RegvTit, group=sampleInfo_RegvTit$sampID)


################################ 
#allsamp = c("10d_RvS","12d_RvS", "R_10dv12d", "S_10dv12d", "10dv12d", "RvS", "PvO")

#allsamp = c("Pas_15dv25d","Reg_15dv25d","Tit_15dv25d","Pas_15dvReg_15d","Pas_15dvTit_15d","Reg_15dvTit_15d","Pas_25dvReg_25d","Pas_25dvTit_25d","Reg_25dvTit_25d","15dv25d","PasvReg","PasvTit","RegvTit")
allsamp = c("PasvReg","PasvTit","RegvTit") #c("PasvTit") #

#all count table:
for (sampID in allsamp) {
dgeFull <- get(paste("dgeFull_",sampID, sep=""))

#Create dir for every sampID
dir.create(file.path(outdir_plot,sampID), showWarnings = TRUE, recursive = FALSE)
outdir_plot_samp = file.path(outdir_plot,sampID, sep = "")

#4.4 Data exploration and quality assessment
##### Exercise 4.7 Extract pseudo-counts (ie log2(K+1))
pseudoCounts <- log2(dgeFull$counts+1)
head(pseudoCounts)

##### Exercise 4.8 Histogram for pseudo-counts (sample Cond.WT.Rep.1)
#hist(pseudoCounts[,"484.10.1"])

##### Exercise 4.9 Boxplot for pseudo-counts
pdf(paste(outdir_plot_samp,sampID,"_Boxplot-for-pseudo-counts.pdf", sep=""))
boxplot(pseudoCounts, col="gray", las=3)
dev.off()

##### Skip error #####
#--------
#dev.off()
#--------
#####Exercise 4.10 MA-plots between WT or Mt samples

#par(mfrow=c(1,2))
## WT1 vs WT2
# A values
#avalues <- (pseudoCounts[,1] + pseudoCounts[,2])/2
# M values
#mvalues <- (pseudoCounts[,1] - pseudoCounts[,2])
#plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="treated")
#abline(h=0, col="red")

## Mt1 vs Mt2
# A values
#avalues <- (pseudoCounts[,4] + pseudoCounts[,5])/2
# M values
#mvalues <- (pseudoCounts[,4] - pseudoCounts[,5])
#plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="control")
#abline(h=0, col="red")

#####Exercise 4.11 MDS for pseudo-counts (using limma package)
pdf(paste(outdir_plot_samp,sampID,"_MDS-for-pseudo-counts.pdf", sep=""))
plotMDS(pseudoCounts, col="blue")
dev.off()
#####Exercise 4.12 heatmap for pseudo-counts
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists

#####
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
pdf(paste(outdir_plot_samp,sampID,"_Heatmap-for-pseudo-counts.pdf", sep=""))
cim(sampleDists, color=cimColor, symkey=FALSE)
dev.off()

#4.5 Differential expression analysis
#Exercise 4.13 remove genes with zero counts for all samples
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
head(dgeFull$counts)

#Exercise 4.14 estimate the normalization factors
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples
head(dgeFull$counts)

eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors
normCounts <- cpm(dgeFull)
pseudoNormCounts <- log2(normCounts + 1)

pdf(paste(outdir_plot_samp,sampID,"_Boxplot-for-pseudoNormCounts.pdf", sep=""))
boxplot(pseudoNormCounts, col="gray", las=3)
dev.off()

pdf(paste(outdir_plot_samp,sampID,"_MDS-for-pseudoNormCounts.pdf", sep=""))
plotMDS(pseudoNormCounts,col="blue")
dev.off()

#Exercise 4.15 estimate common and tagwise dispersion
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull

#Exercise 4.16 perform an exact test for the difference in expression between the two conditions “WT” and “Mt”
dgeTest <- exactTest(dgeFull)
dgeTest

#4.6 Independant filtering
#Exercise 4.17 remove low expressed genes

filtData <- HTSFilter(dgeFull)$filteredData
dgeTestFilt <- exactTest(filtData)
dgeTestFilt

#4.7 Diagnostic plot for multiple testing
#Exercise 4.18 plot an histogram of unadjusted p-values

#hist(dgeTest$table[,"PValue"], breaks=50)

#Exercise 4.19 plot an histogram of unadjusted p-values after filtering
pdf(paste(outdir_plot_samp,sampID,"_Hist-unadjusted_p-values-after-filtering.pdf", sep=""))
hist(dgeTestFilt$table[,"PValue"], breaks=50)
dev.off()

#4.8 Inspecting the results
#Exercise 4.20 extract a summary of the differential expression statistics
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
#head(resNoFilt)
##error occured

resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
#head(resFilt)

#Exporting results to CSV files
write.csv(as.data.frame(resNoFilt),file=paste(outdir_work,sampID,"_results_NoFilt.csv", sep=""))
write.csv(as.data.frame(resFilt),file=paste(outdir_work,sampID,"_results_Filt.csv", sep=""))

# APPLIED ADDITIONAL FILTERS
resFilt_PValue <- resFilt[rownames(resFilt$table)[resFilt$table$PValue<0.01],]
#write.csv(as.data.frame(resFilt_PValue),file=paste(outdir_work,sampID,"_results_Filt_PValue-0.01.csv", sep=""))

resFilt_FDR_logFC <- resFilt[rownames(resFilt$table)[resFilt$table$FDR<0.01 & abs(resFilt$table$logFC)>1.5],]
write.csv(as.data.frame(resFilt_FDR_logFC),file=paste(outdir_work,sampID,"_results_Filt_FDR-0.01-logFC1.5.csv", sep=""))

#4.9 Interpreting the DE analysis results
#Exercise 4.25 create a MA plot with 1% differentially expressed genes
pdf(paste(outdir_plot_samp,sampID,"_MAplot-FDR-0.01-Diff-exp-genes.pdf", sep=""))
plotSmear(dgeTestFilt,de.tags = rownames(resFilt$table)[which(resFilt$table$FDR<0.01)])
dev.off()

#Exercise 4.26 create a Volcano plot
volcanoData <- cbind(resFilt$table$logFC, -log10(resFilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")

#head(volcanoData)
pdf(paste(outdir_plot_samp,sampID,"_Volcano-plot.pdf", sep=""))
plot(volcanoData, pch=19)
dev.off()

#Exercise 4.27 transform the normalized counts in log-counts-per-million
y <- cpm(dgeFull, log=TRUE, prior.count = 1)
head(y)

#Exercise 4.28 select 1% differentially expressed genes and produce a heatmap
selY <- y[rownames(resFilt$table)[resFilt$table$FDR<0.01 & abs(resFilt$table$logFC)>1.5],]
selY_100 <- head(selY, n = 100L)

#pdf(paste(outdir_plot_samp,sampID,"_Heatmap-FDR-0.01-Diff-exp-genes.pdf", sep=""))
pdf(paste(outdir_plot_samp,sampID,"_Heatmap-FDR-0.01-logFC1.5_logCPM_top100.pdf", sep=""))
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]
finalHM <- cim(t(selY_100), color=cimColor, symkey=FALSE)

dev.off()

# Write CPM Values with FDR<0.01 and logFC>1.5 cut-offs
#write.csv(as.data.frame(selY),file=paste(outdir_work,sampID,"_results_Filt_FDR-0.01-logFC1.5_logCPM.csv", sep=""))
#write.csv(as.data.frame(t(finalHM$mat)),file=paste(outdir_work,sampID,"_results_Filt_FDR-0.01-logFC1.5_logCPM_top100.csv", sep=""))
#write.csv(as.data.frame(t(finalHM$mat)),file="/Users/varma/Downloads/PasvsTit_results_Filt_FDR-0.01-logFC1.5_CPM.csv") 

###----
#plot(finalHM$ddc, leaflab="none")
#abline(h=10, lwd=2, col="pink")

#geneClust <- cutree(as.hclust(finalHM$ddc), h=10)
#head(geneClust)

#length(unique(geneClust))
#names(which(geneClust==1))

#########################################################################
}

