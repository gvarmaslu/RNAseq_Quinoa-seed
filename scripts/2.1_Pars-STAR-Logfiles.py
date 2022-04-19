#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# 2.1_Pars-annotation-DEG.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

python2.7 2.1_Pars-STAR-Logfiles.py Pars_STARlogfiles_NG-14833_2.0_Quinoa_map.tsv /home/gala0002/proj/proj_dir/NG-14833_2.0_Align-STAR_Quinoa_Genome/



####------


"""

import sys
import re
import os
import tempfile
import commands
import subprocess
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

	def Search_ReMM(self,outfl1,workdir):
		"""
		Calling Search local DB
		"""
		import glob
		with open(workdir+outfl1,'w') as output:
			headr = str("SampName	Number of input reads	Uniquely mapped reads number	Uniquely mapped reads %	Number of reads mapped to multiple loci	Number of reads mapped to multiple loci %	Number of reads mapped to too many loci	Number of reads mapped to too many loci%	Number of reads unmapped: too short	Number of reads unmapped: too short%	Number of reads unmapped: other	Number of reads unmapped: other%")
			output.write(str(headr)+"\n")
			for name in glob.glob(str(workdir+"*/*.final.out")):
				namesp =  name.split("/")
				filenm = namesp[-1]
				Sampnm = namesp[-2]
				eachln = open(name,'r').readlines()
				eachlnout = str(Sampnm+"\t"+eachln[5].split("|")[1].strip()+"\t"+eachln[8].split("|")[1].strip()+"\t"+eachln[9].split("|")[1].strip()+"\t"+eachln[23].split("|")[1].strip()+"\t"+eachln[24].split("|")[1].strip()+"\t"+eachln[25].split("|")[1].strip()+"\t"+eachln[26].split("|")[1].strip()+"\t"+eachln[30].split("|")[1].strip()+"\t"+eachln[31].split("|")[1].strip()+"\t"+eachln[32].split("|")[1].strip()+"\t"+eachln[33].split("|")[1].strip())
				output.write(str(eachlnout+"\n"))
		print "Done seach for ..."
		return None

clF1 = SearchDB().Search_ReMM(sys.argv[1],sys.argv[2])





