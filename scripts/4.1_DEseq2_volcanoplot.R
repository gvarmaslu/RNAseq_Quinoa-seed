##########################
#  devtools::install_github('kevinblighe/EnhancedVolcano')
#https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

#install.packages("pasilla")
#cat  Express_CDS_Pas-Reg-Tit.gene.counts.matrix | awk '{split($1,a,"_"); print a[4]"_"a[5]"_"a[6]"\t"$0}' | awk '{$2=""; print}' | sed s'/__/Pas_rep1/'g > Express_CDS_Pas-Reg-Tit.gene.counts.matrix_pars 
####
library(DESeq2)
library(edgeR)
library(EnhancedVolcano)
library('pasilla')
####
data = read.table("/Volumes/Mac_HD2/proj_dir/Metadata/Express_counts_allsamp.matrix", header=T, row.names=1, com='')

#### sample 1
col_ordering = c(1,2,3,4,5,6,7,8,9,10,11,12)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("Pas", 6), rep("Reg", 6))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","Pas","Reg")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Pas"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Reg"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="Pas", sampleB="Reg", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$padj),])

#### sample 2
col_ordering = c(1,2,3,4,5,6,13,14,15,16,17,18)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("Pas", 6), rep("Tit", 6))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","Pas","Tit")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Pas"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Tit"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="Pas", sampleB="Tit", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$padj),])

#### sample 3
col_ordering = c(7,8,9,10,11,12,13,14,15,16,17,18)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("Reg", 6), rep("Tit", 6))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","Reg","Tit")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Reg"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Tit"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="Reg", sampleB="Tit", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$padj),])

#### sample 4 - Pas_vs_Reg-Tit
col_ordering = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("Pas", 6), rep("Reg-Tit", 12))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","Pas","Reg-Tit")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Pas"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Reg-Tit"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="Pas", sampleB="Reg-Tit", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$padj),])

#### sample 5 - Reg_vs_Pas-Tit

col_ordering = c(7,8,9,10,11,12,1,2,3,4,5,6,13,14,15,16,17,18)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("AReg", 6), rep("Pas-Tit", 12))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","AReg","Pas-Tit")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "AReg"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Pas-Tit"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="AReg", sampleB="Pas-Tit", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$padj),])

#### sample 6 - Tit_vs_Pas-Reg
col_ordering = c(13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = data.frame(conditions=factor(c(rep("ATit", 6), rep("Pas-Reg", 12))))
rownames(conditions) = colnames(rnaseqMatrix)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","ATit","Pas-Reg")
res = results(dds, contrast)
baseMeanA <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "ATit"])
baseMeanB <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Pas-Reg"])
res = cbind(baseMeanA, baseMeanB, as.data.frame(res))
res = cbind(sampleA="ATit", sampleB="Pas-Reg", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$padj),])


#write.table(res, file='Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Pas_vs_Reg.DESeq2.DE_results', sep='	', quote=FALSE)
#write.table(rnaseqMatrix, file='Express_CDS_Pas-Reg-Tit.gene.counts.matrix.Pas_vs_Reg.DESeq2.count_matrix', sep='	', quote=FALSE)
#source("/data/bioinfo/trinityrnaseq-Trinity-v2.6.5/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")

#### plot
#pdf("Trans.isoform.matrix.Larv_Male_Fem.DESeq2.DE_results.Volcano.pdf")
#plot_MA_and_Volcano(rownames(res), log2(res$baseMean+1), res$log2FoldChange, res$padj)
#dev.off()

#######################
#sampl1
#featureNames <- rownames(res)
#logFoldChange <- res$log2FoldChange
#FDR <- res$padj

#######################
#library(EnhancedVolcano)
#######################
#selectLab = rownames(res)[1:40],
######
#### plot
#pdf("/Volumes/Mac_HD2/proj_dir/NG-14833_4.4_Volcanoplots/Volcano_DESeq2_DE.CDS_Pas-v-Reg.pdf")
#pdf("/Volumes/Mac_HD2/proj_dir/NG-14833_4.4_Volcanoplots/Volcano_DESeq2_DE.CDS_Pas-v-Tit.pdf")
#pdf("/Volumes/Mac_HD2/proj_dir/NG-14833_4.4_Volcanoplots/Volcano_DESeq2_DE.CDS_Reg-v-Tit.pdf")
#
pdf("/Volumes/Mac_HD2/proj_dir/NG-14833_4.4_Volcanoplots/Volcano_DESeq2_DE.CDS_Pas_vs_Reg-Tit_Pva-1e-10_FC-2.pdf")
#pdf("/Volumes/Mac_HD2/proj_dir/NG-14833_4.4_Volcanoplots/Volcano_DESeq2_DE.CDS_Reg_vs_Pas-Tit_Pva-1e-2_FC-1.pdf")
#pdf("/Volumes/Mac_HD2/proj_dir/NG-14833_4.4_Volcanoplots/Volcano_DESeq2_DE.CDS_Tit_vs_Pas-Reg.pdf")

######
#https://www.rdocumentation.org/packages/EnhancedVolcano/versions/1.7.16
######
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                #subtitle = "Pas-vs-Reg",
                #subtitle = "Pas-vs-Tit",
                #subtitle = "Reg-vs-Tit",
                ####
                #subtitle = "Pas_vs_Reg-Tit",
                subtitle = "Reg_vs_Pas-Tit",
                #subtitle = "Tit_vs_Pas-Reg",
                selectLab = head(rownames(res),n = 10L),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 1e-10,
                FCcutoff = 2.0,
                pointSize = 3.0,
                labSize = 3,
                labCol = 'black',
                #labFace = 'bold',
                #boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 15,
                legendIconSize = 4,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

dev.off()
######

# plot_Volcano = function(featureNames, logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20, top_gene_labels_show=20) {
#   plot(logFoldChange, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);
#   text(logFoldChange[1:top_gene_labels_show], (-1*log10(FDR))[1:top_gene_labels_show], labels=featureNames[1:top_gene_labels_show], cex= 0.7, pos=3)
# }
# plot_Volcano(featureNames, logFoldChange, FDR);

####


