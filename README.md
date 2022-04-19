
## Transcriptomics analysis of quinoa seed quality 



>  
* Copy this project to your local directory and access script directory

```
git clone https://github.com/gvarmaslu/RNAseq_Quinoa-seed

cd RNAseq_Quinoa-seed/scripts/

```

> 
NOTE: Before running the workflow one should install all dependent software/packages and set necessary environment paths to the working directly. 


* Example python script for Quality Check (QC)

```
python 1.0_QC.py <proj_dir_path>

```

#### Usage of workflow:

> 1.1
* Run script for RNAseq raw data filter with rRNA and Adapter removal steps

```
bash 1.1_rRNA-and-AdaptRM.sh > 1.1_rRNA-and-AdaptRM.sh.log 2>&1

```

> 1.2
* Run script to merge technical replicates 

```
bash 1.2_Merge-lanes.sh > 1.2_Merge-lanes.sh.log 2>&1

```

> 2.0

* Alignment of RNAseq reads with STAR

```
bash 2.0_Align-STAR.sh > 2.0_Align-STAR.sh.log 2>&1

```

* De novo transcriptome assembly using Trinity

```
bash 2.0_Trinity_de-novo-ass.sh > 2.0_Trinity_de-novo-ass.sh.log 2>&1

--
2.0_Trinity_de-novo_Master-ass.sh
2.1.1_samp-info.tsv
2.1.2_Trinity_Genome-guided.py
2.1.2_Trinity_Genome-guided_Master-ass.sh
2.1_Pars-STAR-Logfiles.py

```


> 3.0

* Quantification analysis

```
bash 3.0_Features-count.sh > 3.0_Features-count.sh.log 2>&1

```

* Script for Parsing feature counts 

```
python2.7 3.1_Features_pars_count.py /home/gala0002/proj/proj_dir/

--
3.0_Features-count.sh
3.1_Features_pars_count.py
3.2.1_create-metadata-table.R
3.2.2_Remove-Columns-in-TSV-file.py
3.2_Features-combine_all.R
```

> 4.0 

* Differential Gene or Transcript Expression (DEG/DET) analysis

```
bash 4.0_DET.sh > 4.0_DET.sh.log 2>&1

--
4.1_DEseq2_volcanoplot.R
4.3_EdgeR_selgenes.R

```

> 5.0

* Annotation and parsing of DEG/DET output files

```
bash 5.0_Pars-annotation-DET_DEseq2_all.sh > 5.0_Pars-annotation-DET_DEseq2_all.log 2>&1

--

5.0_Pars-annotation-DET_DEseq2.py
5.0_Pars-annotation-DET_DEseq2_all.sh
5.1_Pars-annotation-GFF-n-Blast_v2.py

```
> 6.0

* SNP Calling using RNAseq data

```
bash 6.0_SNP-calling.sh > 6.0_SNP-calling.sh.log 2>&1

--
6.0_Pars-annotation-Venn-data.sh
6.0_SNP-calling_KisSplice.sh
6.1_SNP-calling_KissDE.R
6.2_SNP-calling_GATK.sh


```

> 7.0

* Annotation of Differential Gene Expression (DEG) analysis output files 
* Example Usage: 7.0_Pars-annotation-DEG.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

```

python 7.0_Pars-annotation-DEG.py infile infile_anno.tsv <proj dir>

--

7.0_Pars-annotation-DE_trans.py
7.2_Pars-annotation-Venn-data.sh
7.3_Pars-annotation-Venn-data-pars.py
7.4_Pars-SNP-genes.py

```

> 8.0

* Gene Ontology terms and KEGG pathway Enrichment Analysis

```
R CMD 8.0_GO.R

--

8.0_KEGG_Pathway-Enrichment_KEGGREST.R
8.1_KEGG_Pathway-Enrichment_gseKEGG.R
8.2_GO_analysis.R

```


### Citation:

When using this workflow, please cite our publication in Frontiers in Plant Science:

#### Transcriptional regulation of quinoa seed quality: Identification of novel candidate genetic markers for increased protein content
Åsa Grimberg* , Ganapathi Varma Saripella, Ritva Ann-Mari Repo-Carrasco Valencia, Therese Bengtsson, Gabriela Renée Alandia Robles and Anders S. Carlsson
