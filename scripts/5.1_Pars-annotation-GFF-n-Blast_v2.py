#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# 5.0_Pars-annotation-GFF-n-Blast.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

#python 5.0_Pars-annotation-GFF-n-Blast_v2.py Trinity-prot_output.blastx_10 Trinity_TransIDs_Blast-AT-Gene-Anno_v3.tsv /home/gala0002/proj/proj_dir/BLAST_Annotation/
#python 5.0_Pars-annotation-GFF-n-Blast_v2.py Quinoa-ASM168347v1-cds_blastn-max1_AT-000001735.4-cds.tsv Quinoa-ASM168347v1-cds_blastn-max1_AT-000001735.4-cds-Gene-Anno.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_5.0_Annotation/

# BASH script previous step
REF="/home/gala0002/proj/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana"
query="/home/gala0002/proj/proj_dir/REF_Genome/Ref_Chenopodium_Quinoa/GCF_001683475.1_ASM168347v1_cds_from_genomic_clean.fa"
outfl="/home/gala0002/proj/proj_dir/NG-14833_5.0_Annotation/Quinoa-ASM168347v1-cds_blastn-max1_AT-000001735.4-cds.tsv"
NCBIBLAST="/bioinfo/BLAST/ncbi-blast-2.9.0+/bin"
MAGICBLAST="/bioinfo/BLAST/ncbi-magicblast-1.5.0/bin"

#nice -n 5 $NCBIBLAST/blastn -db $REF/GCF_000001735.4_TAIR10.1_cds_from_genomic_clean \
-query $query -num_threads 30 \
-max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 1e-05 -out $outfl

#MAGIC-BASH script previous step
outfl="/home/gala0002/proj/proj_dir/NG-14833_5.0_Annotation/Quinoa-ASM168347v1-cds_magicblast_AT-000001735.4-cds.tsv"
#$MAGICBLAST/magicblast -db $REF/GCF_000001735.4_TAIR10.1_cds_from_genomic_clean -query $query -num_threads 40 -outfmt tabular -out $outfl -reftype transcriptome


"""

import sys
import re
import os
import tempfile
import commands
import subprocess
#import subprocess32
from subprocess import *
from subprocess import call

"""
#Query_ID, Subject_ID, %Identity, Alignment_length, Mismatches, Gap_opens, Q.start, Q.end, S.start, S.end, E-value, Bit-score
#Query_ID	Subject_ID	%Identity	Alignment_length	Mismatches	Gap_opens	Q.start	Q.end	S.start	S.end	E-value	Bit-score
"""


class fileHandler:
	def __init__(self):
		self.data = []
		#print "Calling fileHandler constructor"
	def open_file(self,readfl):
		self.rfile = open(readfl,'r').readlines()
		return self.rfile
	def write_file(self,writefl):
		self.wfile = open(writefl,'w')
		return self.wfile

class SearchDB(fileHandler):

	def __init__(self):
		self.data = []
		from collections import defaultdict
		self.ident_ranges_HMBM = defaultdict(list)

	def Search_ReMM(self,readfl1,outfl1,workdir):
		"""
		Calling Search local DB
		"""
		def srchdb1(CDSName):
			DIR="/Volumes/Seagate/Backup-MAC_HD2/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.4_TAIR10.1_genomic.gff"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(CDSName)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				grepout = cmdFls2.strip()
			except:
				False
			return grepout

		def srchdb2(GName):
			DIR="/Volumes/Seagate/Backup-MAC_HD2/proj_dir/REF_Genome/Ref_Arabidopsis-thaliana/GCF_000001735.4_TAIR10.1_genomic.gff"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w GeneID:'"+str(GName)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				grepout = cmdFls2.strip()
			except:
				False
			return grepout
		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			first_line = f1.readlines()
			#output.write(str("TRINITY_De-novo-Transcript-ID"+"\t"+"AT_Sequence-region-ID"+"\t"+"AT_Transcript-ID"+"\t"+"AT_Gene-Annotation"+"\n"))
			output.write(str("Quinoa-ASM168347v1-CDS-ID"+"\t"+"AT_Transcript-ID"+"\t"+"Percent_Identity	Alignment_length	Mismatches	Gap_opens	Q.start	Q.end	S.start	S.end	E-value	Bit-score"+"\t"+"AT_ChrID"+"\t"+"AT_CDS-Annotation"+"\t"+"AT_Gene-Name"+"\t"+"AT_Gene-Position"+"\t"+"AT_Gene-Annotation"+"\n"))
			for lns in first_line:
				lns_sp0 =  lns.strip().split()[0]
				lns_sp1 =  lns.strip().split()[1]
				lns_sp12 =  "_".join(lns_sp1.split("_")[3:5])
				lns_sp3all =  "\t".join(lns.strip().split()[2:])
				#print lns_sp12
				Anno_out = srchdb1(lns_sp12)
				if len(Anno_out)==0:
					output.write(str(lns_sp0+"\t"+lns_sp1+"\t"+lns_sp3all+"\t"+"."+"\t"+"."+"\n"))
				else:
					Anno_cds = Anno_out.split("\t")
					#print Anno_sp[8].split(";")[2].split("GeneID:")[1].split(",")[0]
					geneidno = Anno_cds[8].split(";")[2].split("GeneID:")[1].split(",")[0]
					Anno_out_gn = srchdb2(geneidno)
					Anno_gene = Anno_out_gn.split("\t")
					Anno_gene_pos = "-".join(Anno_gene[3:5])
					Anno_gene_name = Anno_gene[8].split("Name=")[1].split(";")[0]
					#print Anno_gene_name+"\t"+Anno_gene_pos+"\t"+Anno_gene
					output.write(str(lns_sp0+"\t"+lns_sp1+"\t"+lns_sp3all+"\t"+Anno_cds[0]+"\t"+Anno_cds[8]+"\t"+Anno_gene_name+"\t"+Anno_gene_pos+"\t"+Anno_gene[8]+"\n"))
		print "Done seach for ..."
		return None


clF1 = SearchDB().Search_ReMM(sys.argv[1],sys.argv[2],sys.argv[3])





