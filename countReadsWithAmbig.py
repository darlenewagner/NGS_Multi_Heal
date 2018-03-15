#!/usr/bin/python

import sys
import os.path
import argparse
import re
import warnings

## A command-line script for counting reads containing ambiguous nucleotides
## out of total reads in a fastq-formatted file
## input: filePath/reads.fastq --format simple|integer|percent|verbose
## output: single line to terminal (STDOUT)
### requires prinseq/0.13.x or higher and perl/5.22.x or higher 

## Function: A closure for file extension checking

def ext_check(expected_ext, openner):
	def extension(filename):
		if not filename.lower().endswith(expected_ext):
			raise ValueError()
		return openner(filename)
	return extension

parser = argparse.ArgumentParser(description='This script uses prinseq-lite.pl to find count or proportion of reads with N', usage="countReadsWithAmbig.py filepath/reads.fastq --format [simp(le)?/integer/percent/verbose(default)]")

parser.add_argument("fastq", type=ext_check('.fastq', argparse.FileType('r')))

parser.add_argument("--format", default='verbose', type = lambda s : s.lower(), choices=['simple','integer','percent','verbose','s','simp','i','p','v'])

args = parser.parse_args()

fastqFilePath = args.fastq.name

readCount = os.popen("prinseq-lite.pl -stats_info -fastq {} | tail -1 | perl -ne '$_=~s/stats_info\s+reads\s+//g; print;'".format(fastqFilePath)).read()

intReadCount = int(float(readCount))

#print(intReadCount)

ambigCount = os.popen("prinseq-lite.pl -stats_ns -fastq {} | tail -1 | perl -ne '@F=split(/\t/, $_); print \"$F[2]\";'".format(fastqFilePath)).read()

if(ambigCount):
	intAmbigCount = int(float(ambigCount.rstrip()))
else:
	intAmbigCount = 0

#print(intAmbigCount)

splitFileStr = re.split(pattern='/', string=fastqFilePath)
nameIdx = len(splitFileStr) - 1
fileName = splitFileStr[nameIdx]

if(intReadCount < 1):
	intReadCount = 1
	warnings.warn("No postive value for read counts!")

percentAmbig = float(intAmbigCount/intReadCount)*100

# verbose output: readsFile.fastq         Total Reads: XXXXXXX    Reads with N: YYYY  
if ( args.format == 'verbose' or args.format == 'v' ):
	print(fileName, "\tTotal Reads: ", intReadCount, "\tReads with N: ", intAmbigCount)

# integers-only output: readsFile.fastq    XXXXXXX    YYYY
elif ( args.format == 'integer' or args.format == 'i'):
	print(fileName, "\t", intReadCount, "\t", intAmbigCount)

# percent reads with N: readsFile.fastq    ZZ.ZZZZ %
elif ( args.format == 'percent' or args.format == 'p' ):
	print(fileName, "\t", "%6.4f" % percentAmbig, "%")

# integer N count only: YYYY
elif ( args.format == 'simp' or args.format == 's' or args.format == 'simple'):
	print(intAmbigCount)

args.fastq.close()
