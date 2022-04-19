#!/usr/bin/python


"""
#Script to Pars vcf...


#####------Inputs-------
# ./Sript.py --Workdir
# python2.7 3.1_Features_pars_count.py /home/gala0002/proj/proj_dir/

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
import math

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

	def Search_CADD(self,readfl1):
		"""
		Calling Search 
		"""

		INDIR=readfl1+"NG-14833_3.0_Salmon_Quinoa_Genome/"
		OUTDIR=readfl1+"NG-14833_3.1_Salmon_Quinoa_Genome_pars-count/"
		if not os.path.exists(OUTDIR):
			os.makedirs(OUTDIR)
			
		for filename in os.listdir(INDIR):
			#flout = filename.split("_"))
			#with open(INDIR+filename+"/results.xprs",'r') as f1, open(OUTDIR+filename+"_express.count",'w') as output:
			with open(INDIR+filename+"/salmon/quant.sf",'r') as f1, open(OUTDIR+filename+".count",'w') as output:
				first_line = f1.readline().strip()
				for g in f1:
					g1 = g.strip()
					gg = g1.split("\t")
					#for express
					#gg_count = gg[1]+"\t"+gg[4] 
					#for RSEM
					#gg_count = gg[0]+"\t"+str(math.trunc(float(gg[4])))
					#for Salmon
					#gg_count = gg[0]+"\t"+str(math.trunc(float(gg[4])))
					gg_count = gg[0]+"\t"+str(gg[4])
					output.write(str(gg_count+"\n"))
			print filename
		print "Done parsing...."
		return None

# file1: Input positions SNP
# write1: Output file
clF1 = SearchDB().Search_CADD(sys.argv[1])

