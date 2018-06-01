#!/usr/bin/python

import sys
import os
import os.path
import argparse
import re
import logging
import warnings

## A command-line script for counting reads containing ambiguous nucleotides
## out of total reads in a fastq-formatted file
## input: filePath/reads.fastq --format simple|integer|percent|verbose
## output: single line to terminal (STDOUT)
### requires prinseq/0.13.x or higher and perl/5.22.x or higher 

## Function: A closure for file extension checking

def ext_check(expected_ext1, expected_ext2, openner):
	def extension(filename):
		if not (filename.lower().endswith(expected_ext1) or filename.lower().endswith(expected_ext2)):
			raise ValueError()
		return openner(filename)
	return extension

parser = argparse.ArgumentParser(description='Find count or proportion of reads with N', usage="countReadsWithAmbig.py filepath/reads.fastq --format [simp(le)?/integer/percent/verbose(default)]")

parser.add_argument("fastq", type=ext_check('.fastq', 'fastq.gz', argparse.FileType('r')))

parser.add_argument('--internalN', '-inner' , default='N', choices=['Y', 'N'], help="Count internal Ns only (remove 5-prime and 3-prime Ns first)")

parser.add_argument("--format", default='verbose', type = lambda s : s.lower(), choices=['simple','integer','percent','verbose','s','simp','i','p','v'])

args = parser.parse_args()

logger = logging.getLogger("countReadsWithAmbig.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)


fastqFilePath = args.fastq.name

## Simple file folder path extractor function, when format == 0

def getIsolateStr(filePathString, format):
	splitStr = re.split(pattern='/', string=filePathString)
	fileNameIdx = len(splitStr) - 1
	fileString = splitStr[fileNameIdx]
	if(format == 1):
		isolateString = re.split(pattern='\.', string=splitStr[fileNameIdx])
		#print(isolateString)
		if(re.search(pattern='R1_001_prinseq_R1_001', string=isolateString[0]) and (args.rm_ambig == 'N')):
			## If prior prinseq output and remove_ambig == N, simplify output folder
			fileString = re.sub(r'R1_001_prinseq_R1_001', 'minLength', isolateString[0])
		elif(re.search(pattern='_R1_001', string=isolateString[0])):
			fileString = re.sub(r'R1_001', 'prinseq', isolateString[0])
		else:
			fileString = isolateString[0] + '_prinseq'
	elif(format == 0):  ## convert fileString to directory path
		#fileNameIdx = fileNameIdx - 1
		ii = 1
		fileString = splitStr[0]
		while ii < fileNameIdx:
			fileString = fileString + "/" + splitStr[ii]
			ii = ii + 1
	#print(fileString)
	return fileString


myFilePath = getIsolateStr(fastqFilePath, 0)

tempFilePath = myFilePath + "/tempPrinseq/"
goodFilePath = tempFilePath + "/goodTrEndN"
badFilePath = tempFilePath + "/badTrEndN"
runLog = myFilePath + "/prinseq.runtime.log"

if(args.internalN == 'Y'):
	#if(os.path.exists(myFilePath)):
	os.mkdir(tempFilePath)

readCount = None
fullAmbigCount = None
innerAmbigCount = None
gunzFastqFilePath = None

if(re.search(r'\.fastq$', fastqFilePath, flags=re.IGNORECASE)):
	readCount = os.popen("prinseq-lite.pl -stats_info -fastq {} | tail -1 | perl -ne '$_=~s/stats_info\s+reads\s+//g; print;'".format(fastqFilePath)).read()
	fullAmbigCount = os.popen("prinseq-lite.pl -stats_ns -fastq {} | tail -1 | perl -ne '@F=split(/\t/, $_); print \"$F[2]\";'".format(fastqFilePath)).read()
	if(args.internalN == 'Y'):
		os.system("prinseq-lite.pl -fastq {} -out_good {} -out_bad {} -trim_ns_left 1 -trim_ns_right 1 ns_max_p 93 2> {}".format(fastqFilePath, goodFilePath, badFilePath, runLog)) 
		goodFilePath = goodFilePath + "\.fastq"
		innerAmbigCount = os.popen("prinseq-lite.pl -stats_ns -fastq {} | tail -1 | perl -ne '@F=split(/\t/, $_); print \"$F[2]\";'".format(goodFilePath)).read()
elif(re.search(r'fastq\.gz$', fastqFilePath, flags=re.IGNORECASE)):
	os.system("gunzip -c {} > {}".format(fastqFilePath, gunzFastqFilePath))	
	readCount = os.popen("prinseq-lite.pl -stats_info -fastq {} | tail -1 | perl -ne '$_=~s/stats_info\s+reads\s+//g; print;'".format(gunzFastqFilePath)).read()
	fullAmbigCount = os.popen("prinseq-lite.pl -stats_ns -fastq {} | tail -1 | perl -ne '@F=split(/\t/, $_); print \"$F[2]\";'".format(gunzFastqFilePath)).read()
	if(args.internalN == 'Y'):
		os.system("prinseq-lite.pl -fastq {} -out_good {} -out_bad {} -trim_ns_left 1 -trim_ns_right 1 ns_max_p 93 2> {}".format(fastqFilePath, goodFilePath, badFilePath, runLog)) 
		goodFilePath = goodFilePath + "\.fastq"
		innerAmbigCount = os.popen("prinseq-lite.pl -stats_ns -fastq {} | tail -1 | perl -ne '@F=split(/\t/, $_); print \"$F[2]\";'".format(goodFilePath)).read()
	os.system("rm {}".format(gunzFastqFilePath))
else:
	logger.warn("Input file unrecognized format!")
	sys.exit(1)	

os.system("rm -r {}".format(tempFilePath))

intReadCount = int(float(readCount))

ifullAmbigCount = 0

if(fullAmbigCount):
	ifullAmbigCount = int(float(fullAmbigCount.rstrip()))
else:
	ifullAmbigCount = 0

intInAmbigCount = 0

if(innerAmbigCount):
	intInAmbigCount = int(float(innerAmbigCount.rstrip()))
else:
	intInAmbigCount = 0

#print(ifullAmbigCount)

splitFileStr = re.split(pattern='/', string=fastqFilePath)
nameIdx = len(splitFileStr) - 1
fileName = splitFileStr[nameIdx]

if(intReadCount < 1):
	intReadCount = 1
	warnings.warn("No postive value for read counts!")

percentAmbig0 = float(ifullAmbigCount/intReadCount)*100
percentAmbig1 = float(intInAmbigCount/intReadCount)*100

# verbose output: readsFile.fastq         Total Reads: XXXXXXX    Reads with N: YYYY  
if ( args.format == 'verbose' or args.format == 'v' ):
	if(args.internalN == 'Y'):
		print(fileName, "\tTotal Reads: ", intReadCount, "\tAll reads with N: ", ifullAmbigCount, "\tReads with internal N: ", intInAmbigCount)
	else:
		print(fileName, "\tTotal Reads: ", intReadCount, "\tAll reads with N: ", ifullAmbigCount)
# integers-only output: readsFile.fastq    XXXXXXX    YYYY
elif ( args.format == 'integer' or args.format == 'i'):
	if(args.internalN == 'Y'):
		print(fileName, "\t", intReadCount, "\t", ifullAmbigCount, "\t", intInAmbigCount)
	else:
		print(fileName, "\t", intReadCount, "\t", ifullAmbigCount)

# percent reads with N: readsFile.fastq    ZZ.ZZZZ %
elif ( args.format == 'percent' or args.format == 'p' ):
	if(args.internalN == 'Y'):
		print(fileName, "\t", "%6.4f" % percentAmbig0, "%", "\t", "%6.4f" % percentAmbig1, "%")
	else:
		print(fileName, "\t", "%6.4f" % percentAmbig0, "%")

# integer N count only: YYYY
elif ( args.format == 'simp' or args.format == 's' or args.format == 'simple'):
	print(ifullAmbigCount)

args.fastq.close()
