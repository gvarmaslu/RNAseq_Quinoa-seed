#!/bin/bash


#Pars DEG

echo "run script "

#./5.0_Pars-annotation-DET_DEseq2_all.sh > 5.0_Pars-annotation-DET_DEseq2_all.log 2>&1


###############################################
#all samples 

workdir1="/home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-Reg-Tit/"
workdir2="/home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/"
workdir3="/home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Reg-vs-Pas-Tit/"
workdir4="/home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Tit-vs-Pas-Reg/"
###


<<COMMENT

#python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Reg.DESeq2.DE_results Pas_vs_Reg.DESeq2.DE_results_anno.tsv ${workdir1}
python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Reg.DESeq2.DE_results.P5e-2_C1.DE.subset Pas_vs_Reg.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv ${workdir1}

#python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Tit.DESeq2.DE_results Pas_vs_Tit.DESeq2.DE_results_anno.tsv ${workdir1}
python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Tit.DESeq2.DE_results.P5e-2_C1.DE.subset Pas_vs_Tit.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv ${workdir1}

#python2.7 5.0_Pars-annotation-DET_DEseq2.py Reg_vs_Tit.DESeq2.DE_results Reg_vs_Tit.DESeq2.DE_results_anno.tsv ${workdir1}
python2.7 5.0_Pars-annotation-DET_DEseq2.py Reg_vs_Tit.DESeq2.DE_results.P5e-2_C1.DE.subset Reg_vs_Tit.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv ${workdir1}



#############################################


#python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Reg-Tit.DESeq2.DE_results Pas_vs_Reg-Tit.DESeq2.DE_results_anno.tsv ${workdir2}
python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C1.DE.subset Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv ${workdir2}

#python2.7 5.0_Pars-annotation-DET_DEseq2.py AReg_vs_Pas-Tit.DESeq2.DE_results AReg_vs_Pas-Tit.DESeq2.DE_results_anno.tsv ${workdir3}
python2.7 5.0_Pars-annotation-DET_DEseq2.py AReg_vs_Pas-Tit.DESeq2.DE_results.P5e-2_C1.DE.subset AReg_vs_Pas-Tit.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv ${workdir3}

#python2.7 5.0_Pars-annotation-DET_DEseq2.py ATit_vs_Pas-Reg.DESeq2.DE_results ATit_vs_Pas-Reg.DESeq2.DE_results_anno.tsv ${workdir4}
python2.7 5.0_Pars-annotation-DET_DEseq2.py ATit_vs_Pas-Reg.DESeq2.DE_results.P5e-2_C1.DE.subset ATit_vs_Pas-Reg.DESeq2.DE_results.P5e-2_C1.DE.subset_anno.tsv ${workdir4}


#############################################
COMMENT

python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Reg.DESeq2.DE_results.P1e-2_C1.DE.subset Pas_vs_Reg.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv ${workdir1}
python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset Pas_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv ${workdir1}
python2.7 5.0_Pars-annotation-DET_DEseq2.py Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv ${workdir1}

#############################################

python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv ${workdir2}
python2.7 5.0_Pars-annotation-DET_DEseq2.py AReg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset AReg_vs_Pas-Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv ${workdir3}
python2.7 5.0_Pars-annotation-DET_DEseq2.py ATit_vs_Pas-Reg.DESeq2.DE_results.P1e-2_C1.DE.subset ATit_vs_Pas-Reg.DESeq2.DE_results.P1e-2_C1.DE.subset_anno.tsv ${workdir4}


###
python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C0.5.DE.subset Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C0.5.DE.subset_anno.tsv ${workdir2}
echo "all script done...."


python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.DE.subset Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/
python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Reg.DESeq2.DE_results.P5e-2_C0.DE.subset Pas_vs_Reg.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg/
python2.7 5.0_Pars-annotation-DET_DEseq2.py Pas_vs_Tit.DESeq2.DE_results.P5e-2_C0.DE.subset Pas_vs_Tit.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Tit/

################
################################################################################
#python 7.3_Pars-annotation-Venn-data-pars.py Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq.tsv Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv /Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/

python 7.3_Pars-annotation-Venn-data-pars.py Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.DE.subset_anno_AT.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/
python 7.3_Pars-annotation-Venn-data-pars.py Pas_vs_Reg.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv Pas_vs_Reg.DESeq2.DE_results.P5e-2_C0.DE.subset_anno_AT.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg/
python 7.3_Pars-annotation-Venn-data-pars.py Pas_vs_Tit.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv Pas_vs_Tit.DESeq2.DE_results.P5e-2_C0.DE.subset_anno_AT.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Tit/

################

python 7.3_Pars-annotation-Venn-data-pars.py Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.5.DE.subset_anno.tsv Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.5.DE.subset_anno_AT.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/
python 7.3_Pars-annotation-Venn-data-pars.py Pas_vs_Reg.DESeq2.DE_results.P5e-2_C0.5.DE.subset_anno.tsv Pas_vs_Reg.DESeq2.DE_results.P5e-2_C0.5.DE.subset_anno_AT.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg/
python 7.3_Pars-annotation-Venn-data-pars.py Pas_vs_Tit.DESeq2.DE_results.P5e-2_C0.5.DE.subset_anno.tsv Pas_vs_Tit.DESeq2.DE_results.P5e-2_C0.5.DE.subset_anno_AT.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Tit/


python 7.3_Pars-annotation-Venn-data-pars.py Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C0.5.DE.subset_anno.tsv Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C0.5.DE.subset_anno_AT.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/

################## March 14th,2022

python 7.3_Pars-annotation-Venn-data-pars-v2.py Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C0.5.DE.subset_anno.tsv Pas_vs_Reg-Tit.DESeq2.DE_results.P1e-2_C0.5.DE.subset_anno_AT-GIs.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/
python 7.3_Pars-annotation-Venn-data-pars-v2.py Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.5.DE.subset_anno.tsv Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.5.DE.subset_anno_AT-GIs.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/
python 7.3_Pars-annotation-Venn-data-pars-v2.py Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.DE.subset_anno.tsv Pas_vs_Reg-Tit.DESeq2.DE_results.P5e-2_C0.DE.subset_anno_AT-GIs.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_Quinoa_Pas-vs-Reg-Tit/


################

echo "Script done...."



