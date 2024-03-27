#!/usr/bin/python

import sys
import os
import argparse
import logging
import re
import csv
import numbers

## A command-line script for trimming read pairs by different offsets at 3' ends.  Enables adaptive trim based upon 
## decreases in read length.  Does not trim adapters, at least not yet.
## input: filePath/forwardReads.fastq filePath/reverseReads.fastq --trimF int --trimR int --qualityStats Y/N
## output: trimmed reads and log files in output folder 'TrimByPython/forwardReads.cleaned.fastq/'
### requires >= FASTX-Toolkit-0.0.13

def readable_dir(prospective_dir):
	if not os.path.isdir(prospective_dir):
    		raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
	if os.access(prospective_dir, os.R_OK):
		if( not prospective_dir.endswith("/") ):
			prospective_dir = prospective_dir + "/"
		return prospective_dir
	else:
		raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

## checking that 3-prime --trimF/trimR parameters are between 0 and 90, inclusive

def bandwidth_type(x):
	xx = int(x)
	if( xx < 0 ):
		raise argparse.ArgumentTypeError("Minimum trim should be 0 bp")
	elif( xx > 90 ):
		raise argparse.ArgumentTypeError("Maximum trim should be 90 bp")
	return xx

logger = logging.getLogger("fastxQualAdaptTrimmer_R1andR2.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description="trimming forward reads by -t1 and reverse reads by -t2", usage="python fastxTrimmer_R1andR2.py inputPath/reads_R1_001.fastq inputPath/reads_R2_001.fastq -t1 X -t2 Y --outDir outputPath")

## Trim from 3-prime
parser.add_argument('--trimF', '-t1', required=True, type=bandwidth_type, help="Trim from 3-prim of R1 (0 - 90 bp).")
parser.add_argument('--trimR', '-t2', required=True, type=bandwidth_type, help="Trim R2 from 3-prim of R2 (0 - 90 bp)")

## use fastx_quality_stats for length, N-count, and PHRED quality
## if --qualityStats == N, --trimF and --trimR greater than 30 will default to 30
parser.add_argument('--qualityStats', '-qc', default='Y', choices=['Y', 'N'], help="Compute quality .tsv before trim (Y/N)")

parser.add_argument('--trim_5prime', default='N', choices=['Y', 'N'], help="Trim 1 to 3 bp from 5-prime of both R1 and R2 reads.")
parser.add_argument('--firstPos', type=int, default=1, choices=range(0,3), help="Number of 5-prime positions to trim (0 - 3 bp).")

parser.add_argument('forward', type=argparse.FileType('r'))
parser.add_argument('reverse', type=argparse.FileType('r'))

## output folder
parser.add_argument('--outDir', '-D', type=readable_dir, required=True, action='store')

## cleanup of output folder
parser.add_argument('--clean_output', '-clean', default='Y', choices=['Y','N'], help="Delete fail and singleton files: Y or N (use 'N' when additional processing downstream needed)")

## force overwrite of previous output folder
parser.add_argument('--force', default='N', choices=['Y','N'], help="Overwrite previous output from same filename: Y or N")


args = parser.parse_args()

## output folder

outFolder = args.outDir

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

logger.info("Parameters loaded.")

## outFolder = 'TrimByPython/'

intTrimFwd = int(args.trimF)
intTrimRev = int(args.trimR)


forward = args.forward.name
reverse = args.reverse.name

## two pipelines, one for fastq input, the other for fastq.gz

origWD = os.getcwd()

## Simple filename extractor function

def getIsolateStr(filePathString):
	splitStr = re.split(pattern='/', string=filePathString)
	fileNameIdx = len(splitStr) - 1
	isolateString = re.split(pattern='\.', string=splitStr[fileNameIdx])
	if(re.search(pattern='_R1_001', string=isolateString[0])):
		isolateString = re.sub(r'_R1_001', '_fastxTrim', isolateString[0])
	elif(re.search(pattern='_R1_001', string=isolateString[0]) and re.search(pattern='clean', string=isolateString[1])):
		isolateString = re.sub(r'_R1_001', '_fastxTrim', isolateString[0])
	else:
		isolateString = isolateString[0] + '_fastxTrim'

	if(re.search(pattern='_R1_001_prinseq', string=isolateString)):
		isolateString = re.sub(r'_R1_001_prinseq', '_prinseq', isolateString)
#	if(re.search(pattern='prinseq_R1_001', string=isolateString)):
#		isolateString = re.sub(r'prinseq_R1_001', 'prinseq', isolateString)

	return isolateString

## Simple output file renaming

outputFileString = getIsolateStr(forward)

newOutputFolder = outFolder + outputFileString

## Test if output folder exists
if(os.path.exists(newOutputFolder)):
	if(args.force == 'Y'):
		os.system("rm -v {}/*".format(newOutputFolder))
		os.rmdir(newOutputFolder)
		os.mkdir(newOutputFolder)
	else:
		print("Output folder {} already exists,\n for outDir = {} . . . ".format(outputFileString, outFolder))
		print("Use --force Y or change path for --outDir\nExiting.")
		sys.exit(1)
else:
	os.mkdir(newOutputFolder)


logger.info("Output folder validated/created.")

## requires csv module
def findAvgQual(inFile):
        qualFile = open(inFile, 'r')
        fileText = csv.reader(qualFile, delimiter='\t')
        tableText = [row for row in fileText]
        #print(tableText[0])
        line = 1
        readCount = 0
        minRead = 0
        while line < len(tableText):
                currReadCount = int(tableText[line][1])
                if(currReadCount < readCount):
			#print(tableText[line][0] + " " + tableText[line][1])
                        minRead = int(tableText[line][0]) - 1
                        break
                readCount = currReadCount
                line = line + 1
        #print(currReadCount, " ", readCount, " ", minRead)
        qualFile.close()
        return(minRead)

gunzipForward = re.sub(r'\.gz$', '', forward)
gunzipReverse = re.sub(r'\.gz$', '', reverse)

qualForward = newOutputFolder + "/" + outputFileString + "_R1.quality.tsv"
qualReverse = newOutputFolder + "/" + outputFileString + "_R2.quality.tsv"
fwdMinRead = 0
revMinRead = 0

if(re.search(r'\.fastq$', forward, flags=re.IGNORECASE) and re.search(r'\.fastq$', reverse, flags=re.IGNORECASE)):
	#print("Nice files for", newOutputFolder)
	outputForward = newOutputFolder + "/" + outputFileString + "_R1_001.cleaned.fastq.gz"
	outputReverse = newOutputFolder + "/" + outputFileString + "_R2_001.cleaned.fastq.gz"
	
	if(args.qualityStats == 'Y'):
                logger.info("Beginning fastx_quality_stats")
                os.system("fastx_quality_stats -Q33 -i {} -o {}".format(forward, qualForward))
                os.system("fastx_quality_stats -Q33 -i {} -o {}".format(reverse, qualReverse))
                #os.system("head -2 {}".format(qualForward))
                fwdMinRead = findAvgQual(qualForward)
                revMinRead = findAvgQual(qualReverse)
                if(intTrimFwd < (fwdMinRead - 5)):
                        intTrimFwd = fwdMinRead
                        #print(intTrimFwd)
                if(intTrimRev < (revMinRead - 5)):
                        intTrimRev = revMinRead
                        #print(intTrimRev)
	else:
		if(intTrimFwd > 30):
			intTrimFwd = 30
		if(intTrimRev > 30):
			intTrimRev = 30

	if(args.trim_5prime == 'Y'):
		int5prime = int(args.firstPos) + 1
		logger.info("Beginning 5-prime trim")
		interForward = newOutputFolder + "/prT_" + outputFileString + "_1.fastq" 
		interReverse = newOutputFolder + "/prT_" + outputFileString + "_2.fastq"
		os.system("fastx_trimmer -Q33 -f {} -i {} -o {}".format(int5prime, forward, interForward))
		logger.info("forward (R1) fastq 5-prime trimming complete")
		os.system("fastx_trimmer -Q33 -f {} -i {} -o {}".format(int5prime, reverse, interReverse))
		logger.info("reverse (R2) fastq 5-prime trimming complete")
		os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimFwd, interForward, outputForward))
		logger.info("forward (R1) fastq 3-prime trim complete")
		os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimRev, interReverse, outputReverse))
		logger.info("reverse (R2) fastq 3-prime trim complete")
		os.system("rm -v {}".format(interForward))
		os.system("rm -v {}".format(interReverse))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'fastx_trimmer: All reads trimmed {} bp at 5-prime; Forward reads trimmed {} bp at 3-prime and reverse reads trimmed {} bp at 3-prime' > {}".format(int5prime, intTrimFwd, intTrimRev, outputLog))
	else:	
		os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimFwd, forward, outputForward))
		logger.info("forward (R1) fastq 3-prime trimming complete")
		os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimRev, reverse, outputReverse))
		logger.info("reverse (R2) fastq 3-prime trimming complete")
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'fastx_trimmer: Forward reads trimmed {} bp at 3-prime and reverse reads trimmed {} bp at 3-prime' > {}".format(intTrimFwd, intTrimRev, outputLog))

elif(re.search('\.fastq.gz', forward, flags=re.IGNORECASE) and re.search('\.fastq.gz', reverse, flags=re.IGNORECASE)):
	outputForward = newOutputFolder + "/" + outputFileString + "_R1_001.cleaned.fastq.gz"
	outputReverse = newOutputFolder + "/" + outputFileString + "_R2_001.cleaned.fastq.gz"

	logger.info("gunzip {}".format(forward))
	os.system("gunzip -c {} > {}".format(forward, gunzipForward))
	logger.info("forward (R1) gunzip complete")

	logger.info("gunzip {}".format(reverse))
	os.system("gunzip -c {} > {}".format(reverse, gunzipReverse))
	logger.info("reverse (R2) gunzip complete; start trimming")


	if(args.qualityStats == 'Y'):
		logger.info("Beginning fastx_quality_stats")
		os.system("fastx_quality_stats -Q33 -i {} -o {}".format(gunzipForward, qualForward))
		os.system("fastx_quality_stats -Q33 -i {} -o {}".format(gunzipReverse, qualReverse))
		fwdMinRead = findAvgQual(qualForward)
		revMinRead = findAvgQual(qualReverse)
		if(intTrimFwd > (fwdMinRead - 5)):
			intTrimFwd = fwdMinRead - 5
		if(intTrimRev > (revMinRead - 5)):
			intTrimRev = revMinRead - 5
	else:
		if(intTrimFwd > 30):
			intTrimFwd = 30
		if(intTrimRev > 30):
			intTrimRev = 30


	if(args.trim_5prime == 'Y'):
		int5prime = int(args.firstPos) + 1
		interForward = newOutputFolder + "/prT_" + outputFileString + "_1.fastq"
		os.system("fastx_trimmer -Q33 -f {} -i {} -o {}".format(int5prime, gunzipForward, interForward))
		logger.info("forward (R1) fastq 5-prime trimming complete")
		os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimFwd, interForward, outputForward))
		logger.info("forward (R1) fastq 3-prime trimming complete")
		interReverse = newOutputFolder + "/prT_" + outputFileString + "_2.fastq"
		os.system("fastx_trimmer -Q33 -f {} -i {} -o {}".format(int5prime, gunzipReverse, interReverse))
		logger.info("reverse (R2) fastq 5-prime trimming complete")
		os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimRev, interReverse, outputReverse))
		logger.info("reverse (R2) fastq 3-prime trimming complete")
		logger.info("reverse (R2) fastq file trim complete")
		os.system("rm -v {}".format(interForward))
		os.system("rm -v {}".format(interReverse))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'fastx_trimmer: All reads trimmed {} bp at 5-prime; Forward reads trimmed {} bp at 3-prime and reverse reads trimmed {} bp at 3-prime' > {}".format(int5prime, intTrimFwd, intTrimRev, outputLog))
		logger.info("reverse (R2) fastq file trimming complete")
	else:
		#logger.info("gunzip {}".format(forward))
		#os.system("gunzip -c {} > {}".format(forward, gunzipForward))
		#logger.info("forward (R1) gunzip complete; start trimming")
		os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimFwd, gunzipForward, outputForward))
		#logger.info("forward (R1) fastq file trim complete")
		#logger.info("gunzip {}".format(reverse))
		#os.system("gunzip -c {} > {}".format(reverse, gunzipReverse))
		#logger.info("reverse (R2) gunzip complete; start trimming")
		os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimRev, gunzipReverse, outputReverse))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'fastx_trimmer: Forward reads trimmed {} bp at 3-prime and reverse reads trimmed {} bp at 3-prime' > {}".format(intTrimFwd, intTrimRev, outputLog))
		logger.info("reverse (R2) fastq file trim complete")

	os.system("rm -v {}".format(gunzipForward))
	os.system("rm -v {}".format(gunzipReverse))
else:
	logger.warn("Filetype - R1 {} R2 {} must both end in .fastq".format(forward, reverse))
	logger.warn("Filetype - R1 {} and R2 {} must both end in .fastq.gz".format(forward, reverse))
	logger.exception("Input files must be FASTQ ending in either .fastq or .fastq.gz")


		





