#!/usr/bin/python

import sys
import os
import argparse
from argparse import ArgumentTypeError, ArgumentParser
from pathlib import Path
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

## Check names of input files for .fastq(.gz) or .fq(.gz) suffix
def ext_check(expected_ext1, expected_ext2, expected_ext3, expected_ext4, openner):
	def extension(filename):
		if not (filename.lower().endswith(expected_ext1) or \
                        filename.lower().endswith(expected_ext2) or \
                        filename.lower().endswith(expected_ext3) or \
                        filename.lower().endswith(expected_ext4)):
			raise ValueError()
		return openner(filename)
	return extension
        
## check_for_empty throws an exception when attempting to read empty files
def check_for_empty(fastq_file):
        path1 = Path(fastq_file)
        if(path1.stat().st_size == 0):
            raise ArgumentTypeError(f'{fastq_file} cannot be empty')

## checking that 3-prime --trimF/trimR parameters are between 0 and 25, inclusive
def bandwidth_type(x):
	xx = int(x)
	if( xx < 0 ):
		raise argparse.ArgumentTypeError("Minimum trim should be 0 bp")
	elif( xx > 30 ):
		raise argparse.ArgumentTypeError("Maximum trim should be 30")
	return xx

logger = logging.getLogger("fastxTrimmer_R1andR2.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description="trimming forward reads by -t1 and reverse reads by -t2", usage="python fastxTrimmer_R1andR2.py inputPath/reads_R1_001.fastq inputPath/reads_R2_001.fastq -t1 X -t2 Y --outDir outputPath")

## Trim from 3-prime
parser.add_argument('--trimF', '-t1', required=True, type=bandwidth_type, help="Trim 0 to 30 bp from 3-prime of R1 reads.")
parser.add_argument('--trimR', '-t2', required=True, type=bandwidth_type, help="Trim 0 to 30 bp from R2 reads.")

parser.add_argument('--trim_5prime', default='N', choices=['Y', 'N'], help="Trim 1 to 3 bp from 5-from of both R1 and R2 reads.")
parser.add_argument('--firstPos', type=int, default=1, choices=range(1,4), help="Number of 5-prime positions to trim.")

parser.add_argument('forward', type=ext_check('.fastq', 'fastq.gz', 'fq', 'fq.gz', argparse.FileType('r')))
parser.add_argument('reverse', type=ext_check('.fastq', 'fastq.gz', 'fq', 'fq.gz', argparse.FileType('r')))

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

intTrimFwd = args.trimF
intTrimRev = args.trimR

forward = args.forward.name
reverse = args.reverse.name

## Throw exception for zero-length files
check_for_empty(args.forward.name)
check_for_empty(args.reverse.name)

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

gunzipForward = re.sub(r'\.gz$', '', forward)
gunzipReverse = re.sub(r'\.gz$', '', reverse)

if(re.search(r'\.f(ast)?q$', forward, flags=re.IGNORECASE) and re.search(r'\.f(ast)?q$', reverse, flags=re.IGNORECASE)):
	#print("Nice files for", newOutputFolder)
	outputForward = newOutputFolder + "/" + outputFileString + "_R1_001.cleaned.fastq.gz"
	outputReverse = newOutputFolder + "/" + outputFileString + "_R2_001.cleaned.fastq.gz"

	if(args.trim_5prime == 'Y'):
		int5prime = int(args.firstPos) + 1
		logger.info("Beginning 5-prime trim")
		interForward = newOutputFolder + "/prT_" + outputFileString + "_1.fastq" 
		interReverse = newOutputFolder + "/prT_" + outputFileString + "_2.fastq"
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -f {} -i {} -o {}".format(int5prime, forward, interForward))
		logger.info("forward (R1) fastq 5-prime trimming complete")
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -f {} -i {} -o {}".format(int5prime, reverse, interReverse))
		logger.info("reverse (R2) fastq 5-prime trimming complete")
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimFwd, interForward, outputForward))
		logger.info("forward (R1) fastq 3-prime trim complete")
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimRev, interReverse, outputReverse))
		logger.info("reverse (R2) fastq 3-prime trim complete")
		os.system("rm -v {}".format(interForward))
		os.system("rm -v {}".format(interReverse))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'singularity exec my_fastx_toolkit.sif fastx_trimmer: All reads trimmed {} bp at 5-prime; Forward reads trimmed {} bp at 3-prime and reverse reads trimmed {} bp at 3-prime' > {}".format(int5prime, intTrimFwd, intTrimRev, outputLog))
	else:	
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimFwd, forward, outputForward))
		logger.info("forward (R1) fastq 3-prime trimming complete")
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimRev, reverse, outputReverse))
		logger.info("reverse (R2) fastq 3-prime trimming complete")
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'singularity exec my_fastx_toolkit.sif fastx_trimmer: Forward reads trimmed {} bp at 3-prime and reverse reads trimmed {} bp at 3-prime' > {}".format(intTrimFwd, intTrimRev, outputLog))

elif(re.search(r'\.f(ast)?q.gz', forward, flags=re.IGNORECASE) and re.search(r'\.f(ast)?q.gz', reverse, flags=re.IGNORECASE)):
	outputForward = newOutputFolder + "/" + outputFileString + "_R1_001.cleaned.fastq.gz"
	outputReverse = newOutputFolder + "/" + outputFileString + "_R2_001.cleaned.fastq.gz"

	if(args.trim_5prime == 'Y'):
		int5prime = int(args.firstPos) + 1
		logger.info("gunzip {}".format(forward))
		os.system("gunzip -c {} > {}".format(forward, gunzipForward))
		logger.info("forward (R1) gunzip complete; start trimming")
		interForward = newOutputFolder + "/prT_" + outputFileString + "_1.fastq"
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -f {} -i {} -o {}".format(int5prime, gunzipForward, interForward))
		logger.info("forward (R1) fastq 5-prime trimming complete")
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimFwd, interForward, outputForward))
		logger.info("forward (R1) fastq 3-prime trimming complete")
		logger.info("forward (R1) fastq file trim complete")
		logger.info("gunzip {}".format(reverse))
		os.system("gunzip -c {} > {}".format(reverse, gunzipReverse))
		logger.info("reverse (R2) gunzip complete; start trimming")
		interReverse = newOutputFolder + "/prT_" + outputFileString + "_2.fastq"
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -f {} -i {} -o {}".format(int5prime, gunzipReverse, interReverse))
		logger.info("reverse (R2) fastq 5-prime trimming complete")
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimRev, interReverse, outputReverse))
		logger.info("reverse (R2) fastq 3-prime trimming complete")
		logger.info("reverse (R2) fastq file trim complete")
		os.system("rm -v {}".format(interForward))
		os.system("rm -v {}".format(interReverse))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'singularity exec my_fastx_toolkit.sif fastx_trimmer: All reads trimmed {} bp at 5-prime; Forward reads trimmed {} bp at 3-prime and reverse reads trimmed {} bp at 3-prime' > {}".format(int5prime, intTrimFwd, intTrimRev, outputLog))
		logger.info("reverse (R2) fastq file trimming complete")
	else:
		logger.info("gunzip {}".format(forward))
		os.system("gunzip -c {} > {}".format(forward, gunzipForward))
		logger.info("forward (R1) gunzip complete; start trimming")
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimFwd, gunzipForward, outputForward))
		logger.info("forward (R1) fastq file trim complete")
		logger.info("gunzip {}".format(reverse))
		os.system("gunzip -c {} > {}".format(reverse, gunzipReverse))
		logger.info("reverse (R2) gunzip complete; start trimming")
		os.system("singularity exec my_fastx_toolkit.sif fastx_trimmer -Q33 -t {} -i {} -z -o {}".format(intTrimRev, gunzipReverse, outputReverse))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'singularity exec my_fastx_toolkit.sif fastx_trimmer: Forward reads trimmed {} bp at 3-prime and reverse reads trimmed {} bp at 3-prime' > {}".format(intTrimFwd, intTrimRev, outputLog))
		logger.info("reverse (R2) fastq file trim complete")

	os.system("rm -v {}".format(gunzipForward))
	os.system("rm -v {}".format(gunzipReverse))
else:
	logger.warning("Filetype - R1 {} R2 {} must both end in .fastq".format(forward, reverse))
	logger.warning("Filetype - R1 {} and R2 {} must both end in .fastq.gz".format(forward, reverse))
	logger.exception("Input files must be FASTQ ending in either .fastq or .fastq.gz")


		





