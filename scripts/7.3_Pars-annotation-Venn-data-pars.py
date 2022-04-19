#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# 7.0_Pars-annotation-DEG.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

#python 7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_RegvTit_145_TIT_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_RegvTit_145_TIT_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/
#python 7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_2143_PAS_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_2143_PAS_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/
#python 7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvReg_Inter_RegvTit_165_REG_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvReg_Inter_RegvTit_165_REG_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/

#python 7.3_Pars-annotation-Venn-data-pars.py venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_Inter_RegvTit_PAS-TIT-REG_anno_uniq.tsv venn_Pas-Reg-Tit.DESeq2.DE_results.P1e-2_C1.5_PasvTit_Inter_PasvReg_Inter_RegvTit_PAS-TIT-REG_anno_uniq_pars.tsv /Volumes/Seagate/Backup-MAC_HD2/proj_dir/NG-14833_6.0_Venn-diagram/

python 7.3_Pars-annotation-Venn-data-pars-v2.py Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq.tsv Reg_vs_Tit.DESeq2.DE_results.P1e-2_C1.DE.subset_anno_uniq_pars.tsv /Volumes/Mac_HD2/proj_dir/NG-14833_7.0_filter/

"""

import sys
import re
import os
import tempfile
#import commands
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
		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			first_line0 = f1.readline().strip().split("\t")
			first_lines = f1.readlines() #[1:]
			output.write(str("\t".join(first_line0)+"\t"+"AT_hit"+"\n"))

			for lns in first_lines:
				TF_fam=[]; AT_GID=[]
				lns_sp =  lns.strip().split("\t")
				#lns_sp1 =  lns_sp[0]
				lns_sp = [x.strip(' ') for x in lns_sp]

				#print lns_sp[46], lns_sp[49]
				if (lns_sp[48] != "." ):
					#Anno_out3 = srchdb3(PIDhit1)
					TF_fam.append(str(lns_sp[48].split("|")[1]))
					#print str(lns_sp[46].split("|")[1]+"_family[PlantTFDB/ITAK]")
				elif (lns_sp[51] != "." ):
					#Anno_out3 = srchdb3(PIDhit1)
					TF_fam.append(str(lns_sp[51].split("|")[1]))
					#print str(lns_sp[49].split("|")[1]+"_family[PlantTFDB/ITAK]")

				
				try:
					GI_Desc = lns_sp[46].split("=")[1]
					GI_Desc_pars = GI_Desc.split(",")
					GID=[]
					for i in GI_Desc_pars:
						if (i.split(":")[0] =="GeneID"):
							GID.append(i.split(":")[1])
					Anno_outall = str("\t".join(lns_sp))+"\t"+str(GID[0])+"\n"
				except:
					Anno_outall = str("\t".join(lns_sp))+"\t"+str("No_AT_hit")+"\n"

				output.write(str(Anno_outall))

				"""
				try:
					GI_Desc_pars = GI_Desc.split("%2C")[0]
					
					if (GI_Desc_pars.split()[0] == "uncharacterized" and lns_sp[47].split("=")[0] != "."):
						GI_Desc = GI_Desc_pars+"(AT:"+lns_sp[47].split("=")[1]+")"
					else:
						GI_Desc = GI_Desc_pars
				except:
					False
				if (lns_sp[46] != "." ):
					True
					try:
						AT_GID.append(lns_sp[46].split("Araport:")[1].split(",")[0])
					except IndexError:
						AT_GID.append(lns_sp[46].split("GeneID:")[1].split(",")[0])
					#print str(AT_GID[0])
					if (TF_fam != [] ):
						True
						#print str(GI_Desc+"("+str(AT_GID[0])+")@"+TF_fam[0]+str("_family[PlantTFDB/ITAK]"))
						Anno_outall = str("\t".join(lns_sp))+"\t"+str(GI_Desc+"("+str(AT_GID[0])+")@"+TF_fam[0]+str("_family[PlantTFDB/ITAK]"))+"\t"+str(AT_GID[0])+"\n"
					else:
						#print str(GI_Desc+"("+str(AT_GID[0])+")")
						Anno_outall = str("\t".join(lns_sp))+"\t"+str(GI_Desc+"("+str(AT_GID[0])+")")+"\t"+str(AT_GID[0])+"\n"
				else:
					#print str("No_AT_hit")
					#print GI_Desc
					if (TF_fam != [] ):
						True
						Anno_outall = str("\t".join(lns_sp))+"\t"+str(GI_Desc)+"@"+TF_fam[0]+str("_family[PlantTFDB/ITAK]")+"\t"+str("No_AT_hit")+"\n"
					else:
						Anno_outall = str("\t".join(lns_sp))+"\t"+str(GI_Desc)+"\t"+str("No_AT_hit")+"\n"
				output.write(str(Anno_outall))
				"""

		print("Done seach for ...")
		return None

clF1 = SearchDB().Search_ReMM(sys.argv[1],sys.argv[2],sys.argv[3])





