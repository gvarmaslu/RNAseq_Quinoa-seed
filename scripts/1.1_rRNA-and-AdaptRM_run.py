#!/usr/bin/python

"""
#Script to Pars vcf...
#####------Inputs-------
# python 1.1_rRNA-and-AdaptRM_run.py /home/proj/
"""

import sys
import re
import fnmatch, os
import tempfile
import commands
import subprocess
#import subprocess32
from subprocess import *
from subprocess import call

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
	def Search_CADD(self,readdir):
		"""
		Calling Search localsearch method
		"""
		SCR="/home/gala0002/proj/RNAseq-analysis-run/scripts_quinoa/"
		INDIR=readdir+"Asa-Grimberg_Quinoa-GATC-2018"
		OUTDIR=readdir+"NG-14833_1.1_sort-trim"
		if not os.path.exists(OUTDIR):
			os.makedirs(OUTDIR)
		#cmdFls1 = SCR+"1.1_rRNA-and-AdaptRM.sh NG-14833_Pas_M_1_lib233714_5767_1"
		#grepout1 =  subprocess.check_output(cmdFls1, shell=True)
		#print grepout1
		file2 = self.write_file(SCR+"1.1_rRNA-and-AdaptRM.sh.log")
		for file in os.listdir(INDIR):
			if fnmatch.fnmatch(file, 'NG-14833_*.fastq.gz'):
				fl_spid = file.split("_")
				fl_samid = "_".join(fl_spid[:-2])
				fl_samlib = "_".join(fl_spid[:-1])
				cmdFls1 = SCR+"1.1_rRNA-and-AdaptRM.sh "+fl_samlib
				grepout1 = subprocess.check_output(cmdFls1, shell=True)
				file2.write(str(grepout1)+"\n")
		print "Done parsing...."
		return None

# file1: Input positions SNP
# write1: Output file

clF1 = SearchDB().Search_CADD(sys.argv[1])



