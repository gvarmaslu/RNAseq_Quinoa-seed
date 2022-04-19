library(edgeR)
library(DESeq2)
library(vsn)
library("pheatmap")

#https://molbiocloud.com/help/tiki-index.php?page=Differential-Expression-Analysis-Using-DESeq2#Differential_Expression_Analysis_Using_DESeq2
#https://shiring.github.io/rna-seq/deseq2/teaching/2016/09/29/DESeq2-course
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-sample-to-sample-distances

#res = read.table("/Volumes/Seagate/Backup-MAC_HD2/proj_Santosh/scripts/Trans.isoform.counts.matrix.AtgamFem_vs_AtgamMale.DESeq2.DE_results",sep="\t", header = TRUE, row.names = 1)
#res[ , !(names(res) %in% c("sampleA","sampleB","baseMeanA","baseMeanB") )]

#Importing the files  
#outdir_plot="/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_4.3_DEG_results/Results_plots/"
outdir_plot="/Volumes/Mac_HD2/proj_dir/NG-14833_8.0_heatmaps/"

dir.create(file.path(outdir_plot), showWarnings = F, recursive = FALSE)

##################################################################################################
#allsamp_UP = c("AtgamFem-UP", "AtgamLarv-UP", "AtgamMale-UP")
#for (sampID_UP in allsamp_UP) {

#allsamp = c("Pas_vs_Reg", "Pas_vs_Tit", "Reg_vs_Tit")
#allsamp = c("Pas_vs_Reg")
#allsamp = c("Pas_vs_Tit")
#allsamp = c("Reg_vs_Tit")

#allsamp = c("Pas_Reg_Tit")
#allsamp = c("REG")
#allsamp = c("PAS")
#allsamp = c("TIT")

allsamp = c("Pas_vs_Reg-Tit")
#allsamp = c("Reg_vs_Pas-Tit")
#allsamp = c("Tit_vs_Pas-Reg")
#all count table:
for (sampID in allsamp) {
  dgeFull_ann <- "/Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C3.DE.subset_anno_uniq_pars_sel.tsv"
  #dgeFull_ann <- "/Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/Reg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv"
  #dgeFull_ann <- "/Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/Tit_vs_Pas-Reg.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv"
  
  #dgeFull_ann <- "/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvReg_1371_anno_uniq_pars.tsv"
  #dgeFull_ann <- "/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_645_anno_uniq_pars.tsv"
  #dgeFull_ann <- "/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_RegvTit_96_anno_uniq_pars.tsv"
  
  #dgeFull_ann <- "/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_Inter_RegvTit_PAS-TIT-REG_anno_uniq_pars.tsv"
  #dgeFull_ann <- "/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvReg_Inter_RegvTit_165_REG_anno_uniq_pars.tsv"
  #dgeFull_ann <- "/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_2143_PAS_anno_uniq_pars.tsv"
  #dgeFull_ann <- "/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_RegvTit_145_TIT_anno_uniq_pars.tsv"
  
  dgeFull_mat <- "/Volumes/Mac_HD2/proj_dir/Metadata/Express_counts_allsamp.matrix"
  #dgeFull_mat <- paste("/Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_2.1_STAR_Genome-CDS-ass/DESeq2_genes/Express_CDS_Pas-Reg-Tit.gene.counts.matrix.",sampID,".DESeq2.count_matrix", sep="")
  
  #res_sel = read.table(dgeFull_ann,sep="\t", header = TRUE, row.names = 1, fill=FALSE)
  res_sel = read.table(dgeFull_ann,sep="\t", header = TRUE, row.names = 1, comment.char="#", na.strings=".", stringsAsFactors=FALSE, quote="", fill=FALSE)
  
  read.mat = read.table(dgeFull_mat,sep="\t", header = TRUE, row.names = 1)
  #coldata = data.frame(row.names = colnames(read.mat),group = rep(c("gt1","gt2"),each = 6),treatment = rep(c(unlist(strsplit(sampID, "_vs_"))),each = 6))
  
  coldata = data.frame(row.names = colnames(read.mat),treatment = c(rep("Pas",each = 6),rep("RegTit",each = 12)) )
  #coldata = data.frame(row.names = colnames(read.mat),treatment = c(rep("PasTit",each = 6),rep("Reg",each = 6),rep("PasTit",each = 6)) )
  #coldata = data.frame(row.names = colnames(read.mat),treatment = c(rep("PasReg",each = 12), rep("Tit",each = 6)) )
  
  ####
  # Define the annotation
  #annotation_row = data.frame(res_sel$Gene.Ann6)
  
  #4. Run DESeqDataSetFromMatrix
  dds <- DESeqDataSetFromMatrix(countData = round(read.mat), colData = coldata, design = ~ treatment) 
  de.genes <- rownames(res_sel)
  
  ####### Select the genes ##### tyoe1
  #dds_sel <- dds[de.genes,]
  #### Stabilizing method 
  #vsd_sel <- normTransform(dds_sel)
  #vsd_sel <- rlog(dds_sel,fitType='mean', blind=FALSE)
  #vsd_sel <- varianceStabilizingTransformation(dds_sel,fitType='local', blind=FALSE)
  #### replace row names
  #vst_sel_ids <- assay(vsd_sel)
  #rownames(vst_sel_ids) <- res_sel$Gene.Ann6
  #mat = vst_sel_ids[ order(res_sel$padj), ] # select the all genes with the lowest padj
  #mat = mat - rowMeans(mat) # Subtract the row means from each value
  
  #######
  #### Stabilizing method 
  #vsd <- normTransform(dds)
  vsd <- rlog(dds, blind=FALSE)
  #vsd <- vst(dds, blind=FALSE)
  #### write.table(assay(vsd), file="/Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C3.DE.subset_anno_uniq_pars-matrix.tsv", sep="\t", quote=FALSE, row.names = T) 
  ######## Select the genes ##### type 2
  vsd_sel <- vsd[de.genes,]
  #### replace row names 
  vsd_sel_ids <- assay(vsd_sel)
  rownames(vsd_sel_ids) <- res_sel$Gene.Ann6
  #### #mat = assay(vsd_sel)[ order(res_sel$padj), ] # select the all genes with the lowest padj
  mat = vsd_sel_ids[ order(res_sel$padj), ] # select the all genes with the lowest padj
  #mat = mat - rowMeans(mat) # Subtract the row means from each value
  
  ###################
  # replace row names 
  #rownames(mat) <- res_sel$Gene.Ann6
  
  # Optional, but to make the plot nicer:
  df = as.data.frame(colData(vsd_sel)[,c("treatment")]) # Create a dataframe with a column of the conditions
  colnames(df) = "condition" # Rename the column header
  rownames(df) = colnames(mat) # add rownames
  ###########
  #scaling
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  data_subset_norm <- t(apply(mat, 1, cal_z_score))
  #pheatmap(data_subset_norm)
  ####################
  ####################
  ## and plot the actual heatmap
  #pheatmap( vsd_sel_ids - rowMeans(vsd_sel_ids) , annotation_col=df)
  #pheatmap(mat, cellwidth = 25, cellheight = 12, fontsize = 12, annotation_col=df)
  pheatmap(mat, cellwidth = 25, cellheight = 12, fontsize = 12, filename = paste(outdir_plot,sampID,"_Heatmap-Diff-exp-genes_selgenes_3cult_P1e-2_C3_rlog.pdf", sep=""), annotation_col=df)
  #########################################################################
}



