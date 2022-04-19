#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# 7.0_Pars-annotation-DET.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

#python 7.4_Pars-SNP-genes.py Tit_vs_Pas-Reg_Genes_uniq_SNPsel_High-impact_Tit.tsv Tit_vs_Pas-Reg_Genes_uniq_SNPsel_High-impact_Tit_pars.tsv //Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/

####------


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
		def srchdb0(GName):
			DIR="/home/gala0002/proj/proj_dir/REF_Genome/Ref_Chenopodium_Quinoa/GCF_001683475.1_ASM168347v1_genomic.gff"
			try:
				True
				#grep -iE "cds.*rna63925" -w -m1 GCF_001683475.1_ASM168347v1_genomic.gff

				cmdFls1 = "LANG=C grep -m1 -w -iE 'cds.*"+str(GName)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout = str(ann[8].split(";protein_id=")[1])
			except:
				False
				grepout = str(".")
			return grepout

		def srchdb1(GName1):
			DIR="/home/gala0002/proj/proj_dir/NG-14833_5.0_Annotation/Quinoa-ASM168347v1-cds-Gene-Anno.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName1)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				grepout1 = str(cmdFls2.strip())
			except:
				False
				grepout1 = str("."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout1

		def srchdb2(GName2):
			DIR="/home/gala0002/proj/proj_dir/NG-14833_5.0_Annotation/working_files/outfile_matchgeneIDs-phmmer-blast_tophits.txt"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName2)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AT_ann = cmdFls2.strip().split()
				#print AT_ann
				grepout2 = str(AT_ann[0]+"\t"+AT_ann[4]+"\t"+AT_ann[5]+"\t"+str(" ".join(AT_ann[18:])))
			except:
				False
				grepout2 = str("."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout2

		def srchdb3(GName3):
			DIR="/home/gala0002/proj/proj_dir/NG-14833_5.0_Annotation/working_files/Quinoa_vs_Thaliana_blastp-max1_proteins_anno.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName3)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AT_ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout3 = str(AT_ann[1]+"\t"+AT_ann[10]+"\t"+AT_ann[2]+"\t"+AT_ann[-2]+"\t"+AT_ann[-1])
			except:
				False
				grepout3 = str("."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout3

		def srchdb4(GName4):
			DIR="/home/gala0002/proj/proj_dir/NG-14833_5.0_Annotation/Quinoa-ASM168347v1-cds_blastn-max1_PlantTFDBv5.0_Ext-Cqu-cds.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 '"+str(GName4)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout4 = str(ann[1]+"\t"+ann[-2]+"\t"+ann[2])
			except:
				False
				grepout4 = str("."+"\t"+"."+"\t"+".")
			return grepout4

		def srchdb5(GName5):
			DIR="/home/gala0002/proj/proj_dir/NG-14833_5.0_Annotation/Quinoa-ASM168347v1-cds_blastn-max1_ITAKv18.12_Ext-Cqu-cds.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 '"+str(GName5)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout5 = str(ann[1]+"\t"+ann[-2]+"\t"+ann[2])
			except:
				False
				grepout5 = str("."+"\t"+"."+"\t"+".")
			return grepout5

		def srchdb6(GName6):
			DIR="/home/gala0002/proj/proj_dir/NG-14833_5.0_Annotation/Quinoa-ASM168347v1-prot_blastp-max1_PI614886-prot.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 '"+str(GName6)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout6 = str(ann[1].split("-")[0]+"\t"+ann[-2]+"\t"+ann[2])
			except:
				False
				grepout6 = str("."+"\t"+"."+"\t"+".")
			return grepout6

		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			#first_line0 = f1.readline().strip().split("\t")
			first_lines = f1.readlines() #[1:]
			head1 = str("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	FORMAT(GT:AD:DP:GQ:PL)	Trans-ID	Up-and-Downstream_Trans-IDs	Gene_ID	Up-and-Downstream_Gene-IDs")
			head2 = str("Quinoa-ASM168347v1-CDS-ID	Protein-ID	Gene_ID	Gene_Name	Gene_Pos	Gene-Ann1	Gene-Ann2	Gene-Ann3	Gene-Ann4	Gene-Ann5	Target_AT_PID_HMMER	Target_AT_E-value_HMMER	Target_AT_Score_HMMER	Target_AT_Description_HMMER	Target_AT_PID_NCBIBLAST	Target_AT_E-value_NCBIBLAST	Target_AT_Score_NCBIBLAST	Target_AT_GeneID_NCBIBLAST	Target_AT_Description_NCBIBLAST	PlantTFDBv5.0_Ext-Cqu-cds_TFfamID	PlantTFDBv5.0_Evalue	PlantTFDBv5.0_Identity_Score	ITAKv18.12_Ext-Cqu-cds_TFfamID	ITAKv18.12_Evalue	ITAKv18.12_Identity_Score	Target_Cq.PI614886_PID_NCBIBLAST	Target_Cq.PI614886_E-value_NCBIBLAST	Target_Cq.PI614886_Score_NCBIBLAST")
			output.write(head1+"\t"+head2+"\n")
			for lns in first_lines:
				lns_sp = lns.strip().split("\t")
				lns_sp_trans = lns_sp[7].split("|transcript|rna")
				lns_sp_genes = lns_sp[7].split("|transcript|XM_")
				allrna =[]; allgenes=[]
				for ii in lns_sp_trans[1:]:
					allrna.append(str("rna"+ii.split("|")[0]))
				firstrna = str("rna")+lns_sp_trans[1].split("|")[0]
				#print firstrna,";".join(allrna)
				#print lns_sp_genes[:-1]
				for ii in lns_sp_genes[:-1]:
					allgenes.append(ii.split("|")[-1])
				firstgns = lns_sp_genes[0].split("|")[-1]
				#print firstgns,";".join(allgenes)
				Anno_out0 = srchdb0(firstrna)
				#print Anno_out0 
				#######
				Anno_out1 = srchdb1(Anno_out0)
				Anno_out2 = srchdb2(Anno_out0)
				Anno_out3 = srchdb3(Anno_out0)
				Anno_out4 = srchdb4(Anno_out0)
				Anno_out5 = srchdb5(Anno_out0)
				Anno_out6 = srchdb6(Anno_out0)
				Anno_outall = str(lns.strip()+"\t"+firstrna+"\t"+";".join(allrna)+"\t"+firstgns+"\t"+";".join(allgenes)+"\t"+Anno_out1+"\t"+Anno_out2+"\t"+Anno_out3+"\t"+Anno_out4+"\t"+Anno_out5+"\t"+Anno_out6+"\n")
				#print Anno_outall
				output.write(str(Anno_outall))
		print "Done seach for ..."
		return None

clF1 = SearchDB().Search_ReMM(sys.argv[1],sys.argv[2],sys.argv[3])





