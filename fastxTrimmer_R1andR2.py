#!/usr/bin/python

import sys
import os
import argparse
import logging
import re

## A command-line script for trimming read pairs by different offsets at 3' ends
## Not dependent upon Bpipe, Java, or Perl
## input: filePath/forwardReads.fastq filePath/reverseReads.fastq --trimF int --trimR int
## output: trimmed reads and log files in output folder 'TrimByPython/forwardReads.cleaned.fastq/'
### requires >= FASTX-Toolkit-0.0.13, >= prinseq-lite.pl 0.20.x, and >= perl/5.22.x

def readable_dir(prospective_dir):
	if not os.path.isdir(prospective_dir):
    		raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
	if os.access(prospective_dir, os.R_OK):
		if( not prospective_dir.endswith("/") ):
			prospective_dir = prospective_dir + "/"
		return prospective_dir
	else:
		raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))


logger = logging.getLogger("fastxTrimmer_R1andR2.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description="trimming forward reads by -t1 and reverse reads by -t2", usage="python fastxTrimmer_R1andR2.py inputPath/reads_R1_001.fastq inputPath/reads_R2_001.fastq -t1 X -t2 Y --outDir outputPath")

parser.add_argument('--trimF', '-t1', required=True, type=int, choices=range(0,21))
parser.add_argument('--trimR', '-t2', required=True, type=int, choices=range(0,21))

parser.add_argument('forward', type=argparse.FileType('r'))
parser.add_argument('reverse', type=argparse.FileType('r'))

parser.add_argument('--outDir', '-D', type=readable_dir, required=True, action='store')

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

intTrimFwd = args.trimF
intTrimRev = args.trimR

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
		isolateString = re.sub(r'R1_001', 'fastxTrim_fastq', isolateString[0])
	else:
		isolateString = isolateString + '_fastxTrim_fastq'
	return isolateString

## Simple output file renaming

outputFileString = getIsolateStr(forward)


newOutputFolder = outFolder + outputFileString

incre = 1
if(os.path.exists(newOutputFolder)):
	while(os.path.exists(newOutputFolder)):
		incre = incre + 1
		newOutputFolder = newOutputFolder + "_" + str(incre)
		if(os.path.exists(newOutputFolder) is False):
			os.mkdir(newOutputFolder)
			break
else:
	os.mkdir(newOutputFolder)
incre = 0

logger.info("Output folder validated/created.")

if(re.search('\.fastq', forward, flags=re.IGNORECASE) and re.search('\.fastq', reverse, flags=re.IGNORECASE)):
	#print("Nice files for", newOutputFolder)
	outputForward = newOutputFolder + "/" + outputFileString + "_R1_001.cleaned.fastq"
	os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimFwd, forward, outputForward))
	logger.info("forward (R1) fastq file trimming complete")
	outputReverse = newOutputFolder + "/" + outputFileString + "_R2_001.cleaned.fastq"
	os.system("fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimRev, forward, outputReverse))
	logger.info("reverse (R2) fastq file trimming complete")
	outputLog = newOutputFolder + "/parameters.log"
	os.system("echo 'fastx_trimmer: Forward reads trimmed {} and reverse reads trimmed {}' > {}".format(intTrimFwd, intTrimRev, outputLog))

 ##elif(re.search('\.fastq.gz', forward, flags=re.IGNORECASE) and re.search('\.fastq.gz', reverse, flags=re.IGNORECASE)):
	





