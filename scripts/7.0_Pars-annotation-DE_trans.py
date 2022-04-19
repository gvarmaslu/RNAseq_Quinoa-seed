#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# 7.0_Pars-annotation-DEG.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

#python 7.0_Pars-annotation-DEG_v1.py Trinity_trans.isoform.counts.matrix.Reg_vs_Tit.DESeq2.DE_results.P1e-3_C2.Tit-UP.subset Trinity_trans.isoform.counts.matrix.Reg_vs_Tit.DESeq2.DE_results.P1e-3_C2.Tit-UP.subset_anno.tsv /home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Master_Tit-Reg-Pas/DESeq2_trans/

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
		def srchdb1(GName1):
			DIR="/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Master_Tit-Reg-Pas/Trinity_blastx-max1-ASM168347v1_protein.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName1)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AT_ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout1 = str("\t".join(AT_ann[1:]))
			except:
				False
				grepout1 = str("."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout1
		def srchdb2(GName2):
			DIR="/home/gala0002/proj/proj_dir/NG-14833_2.0_Trinity_Master_Tit-Reg-Pas/Trinity_blastx-max1-TAIR10.1_protein.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName2)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AT_ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout2 = str("\t".join(AT_ann[1:]))
			except:
				False
				grepout2 = str("."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout2
		def srchdb3(GName3):
			DIR="/home/gala0002/proj/RNAseq-analysis-run/scripts_matchgeneIDs-ncbi-blast/working_files/Quinoa_vs_Thaliana_blastp-max1_proteins_anno.tsv"
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName3)+"' "+str(DIR)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AT_ann = cmdFls2.strip().split("\t")
				#print AT_ann
				grepout3 = str(AT_ann[0]+"-@-"+AT_ann[1]+"\t"+str("\t".join(AT_ann[2:-3]))+"\t"+AT_ann[-3].split()[0]+"\t"+AT_ann[-2]+"\t"+AT_ann[-1])
			except:
				False
				grepout3 = str(GName3+"-@-"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout3

		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			first_line0 = f1.readline().strip().split("\t")
			first_lines = f1.readlines()[1:]
			#output.write(str("Quinoa_GeneID"+"\t"+str("\t".join(first_line0))+"\t"+"Gene_ID	Gene_Name	Gene_Pos	Gene-Ann1	Gene-Ann2	Gene-Ann3	Gene-Ann4	Gene-Ann5	Target_AT_PID_HMMER	Target_AT_E-value_HMMER	Target_AT_Description_HMMER	Target_AT_PID_NCBIBLAST	Target_AT_E-value_NCBIBLAST	Target_AT_GeneID_NCBIBLAST	Target_AT_Description_NCBIBLAST"+"\n"))
			outhdr1 = str("Subject_PID-ASM168347v1_NCBIBLASTX-@-Query_GeneID	Subject_PID-ASM168347v1_NCBIBLASTX	Identity	Alignment_length	Mismatches	Gap_opens	Q.start	Q.end	S.start	S.end	E-value	Bit-score")
			outhdr2 = str("Subject_PID-TAIR10.1_NCBIBLASTX-@-Query_GeneID	Subject_PID-TAIR10.1_NCBIBLASTX	Identity	Alignment_length	Mismatches	Gap_opens	Q.start	Q.end	S.start	S.end	E-value	Bit-score")
			outhdr3 = str("Query_ASM168347v1_GeneID-@-Subject_TAIR10.1_GeneID	Identity	Alignment_length	Mismatches	Gap_opens	Q.start	Q.end	S.start	S.end	E-value	Bit-score	Dbxref	product")
			output.write(str("GeneID"+"\t"+str("\t".join(first_line0))+"\t"+outhdr1+"\t"+outhdr2+"\t"+outhdr3+"\n"))
	
			for lns in first_lines:
				lns_sp =  lns.strip().split("\t")
				lns_sp1 =  lns_sp[0]
				
				#protID = lns_sp1[3]+"_"+lns_sp1[4]
				Anno_out1 = srchdb1(lns_sp1)
				Anno_out2 = srchdb2(lns_sp1)
				Anno_out3 = srchdb3(lns_sp1)

				#print Anno_out1, Anno_out2, len(Anno_out1.split("\t")), len(Anno_out2.split("\t"))
				PIDhit1 = Anno_out1.split("\t")[0]
				PIDhit2 = Anno_out2.split("\t")[0]
				isoform_id_an1 = PIDhit1+"@"+str("_".join(lns_sp1.split("_")[1:]))
				isoform_id_an2 = PIDhit2+"@"+str("_".join(lns_sp1.split("_")[1:]))

				if (PIDhit1 != "."):
					Anno_out3 = srchdb3(PIDhit1) 
				else:
					Anno_out3 = str(PIDhit1+"-@-"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+".")

				Anno_outall = str("\t".join(lns_sp))+"\t"+str(isoform_id_an1)+"\t"+str(Anno_out1)+"\t"+str(isoform_id_an2)+"\t"+str(Anno_out2)+"\t"+str(Anno_out3)+"\n"
				#print str("\t".join(lns_sp))+"\t"+str(Anno_out3)
				output.write(str(Anno_outall))
		print "Done seach for ..."
		return None

clF1 = SearchDB().Search_ReMM(sys.argv[1],sys.argv[2],sys.argv[3])





