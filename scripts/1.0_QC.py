#!/usr/bin/python

"""
#Script to Pars vcf...
#####------Inputs-------
# python 1.0_QC.py /home/proj/
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
		INDIR=readdir+"Asa-Grimberg_Quinoa-GATC-2018"
		OUTDIR=readdir+"NG-14833_1.0_QC"
		if not os.path.exists(OUTDIR):
			os.makedirs(OUTDIR)
		for file in os.listdir(INDIR):
			if fnmatch.fnmatch(file, 'NG-14833_*.fastq.gz'):
				fl_spid = file.split("_")
                                fl_samid = "_".join(fl_spid[:-2])
				fl_samlib = "_".join(fl_spid[:-1])
				print fl_samid, fl_samlib, fl_spid[3]
				cmdFls1 = "cd "+INDIR+"; mkdir -p "+OUTDIR+"/"+fl_samid
				cmdFls2 = "nice -n 5 fastqc -o "+OUTDIR+"/"+fl_samid+" -t 20 --noextract "+fl_samlib+"*_1.fastq.gz "+fl_samlib+"*_2.fastq.gz"
				grepout1 =  subprocess.check_output(cmdFls1+";"+cmdFls2, shell=True)
				print grepout1
		print "Done parsing...."
		return None

# file1: Input positions SNP
# write1: Output file

clF1 = SearchDB().Search_CADD(sys.argv[1])
