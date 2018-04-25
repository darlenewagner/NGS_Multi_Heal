#!/usr/bin/python

import sys
import os
import argparse
import logging
import re

## A command-line script for cleaning read pairs of reads < XX bp and reads with > Y Ns
## Offers length-based and N-presence-based cleaning only, no quality checking

## input: filePath/forwardReads.fastq filePath/reverseReads.fastq --rm_ambig Y/N --outDir (optional: )
## output: trimmed reads and log files in output folder 'TrimByPython/forwardReads_prinseq/'
### requires >= prinseq-lite.pl 0.20.x and >= perl/5.22.x

def readable_dir(prospective_dir):
	if not os.path.isdir(prospective_dir):
    		raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
	if os.access(prospective_dir, os.R_OK):
		if( not prospective_dir.endswith("/") ):
			prospective_dir = prospective_dir + "/"
		return prospective_dir
	else:
		raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

## for checking that --min-length is between 35 and 145, inclusive
def bandwidth_type(x):
	xx = int(x)
	if( xx < 34 ):
		raise argparse.ArgumentTypeError("Minimum min-length should be 35")
	elif( xx > 144 ):
		raise argparse.ArgumentTypeError("Maximum min-length should be 145")
	return xx


logger = logging.getLogger("prinseqLite_R1andR2.py")
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description="Wrapper for prinseq-lite.pl", usage="python simpPrinseqLite_R1andR2.py inputPath/reads_R1_001.fastq inputPath/reads_R2_001.fastq --rm_ambig Y/N --min_len 100  --trim_qual [Y/N] --outDir outputPath")

## minimum read length
parser.add_argument('--min_len', '-min', type=bandwidth_type, default=100, help="Remove reads less than --min_len: Default=100, 39 < --min_len < 146")

## remove Ns
parser.add_argument('--rm_ambig', '-rm', required=True, choices=['Y','N'], help="Remove amp")
parser.add_argument('--ambig_allow', '-ns_max_n', type=int, default=1, choices=range(0,5), help="Ns tolerance: 0 to 5 per read, default = 1")

## Trim by PHRED Q-score
parser.add_argument('--trim_qual', default='N', choices=['Y','N'], help="trim by base pair quality (Y|N)?")
parser.add_argument('--min_score', default=24, type=int, choices=range(19,30), help="quality score on which to trim [20 to 30: default 24]")

## input files
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

intMinLen = args.min_len
intTrimN = args.ambig_allow

forward = args.forward.name
reverse = args.reverse.name

logger.info("Parameters loaded.")

## two pipelines, one for fastq input, the other for fastq.gz

origWD = os.getcwd()

## Simple filename extractor function

def getIsolateStr(filePathString, format):
	splitStr = re.split(pattern='/', string=filePathString)
	fileNameIdx = len(splitStr) - 1
	fileString = splitStr[fileNameIdx]
	if(format == 1):
		isolateString = re.split(pattern='\.', string=splitStr[fileNameIdx])
		#print(isolateString)
		if(re.search(pattern='_R1_001\.fastq', string=isolateString[0])):
			fileString = re.sub(r'R1_001\.fastq', 'prinseq', isolateString[0])
		else:
			fileString = isolateString[0] + '_prinseq'
	#print(fileString)
	return fileString

## Simple output file renaming

outputFileString = getIsolateStr(forward, 1)

newOutputFolder = outFolder + outputFileString

## Validate input file names: FILTER 1

if(forward == reverse):
	logger.warn("Forward reads and reverse reads files cannot have the same name and path!")
	sys.exit(1)
#elif()

firstLineFwd = None
firstLineRev = None
lineCountFwd = 0
lineCountRev = 0

## Validate input file suffixes: FILTER 2

if(re.search(r'\.gz$', forward, flags=re.IGNORECASE) and re.search(r'\.gz$', reverse, flags=re.IGNORECASE)):
	firstLineFwd = os.popen("zcat {} | head -1".format(forward)).read()
	firstLineRev = os.popen("zcat {} | head -1".format(reverse)).read()
	lineCountFwd = os.popen("zcat {} | wc -l".format(forward)).read()
	lineCountRev = os.popen("zcat {} | wc -l".format(reverse)).read()
elif(re.search(r'\.fastq$', forward, flags=re.IGNORECASE) and re.search(r'\.fastq$', reverse, flags=re.IGNORECASE)):
	firstLineFwd = os.popen("cat {} | head -1".format(forward)).read()
	firstLineRev = os.popen("cat {} | head -1".format(reverse)).read()
	lineCountFwd = os.popen("cat {} | wc -l".format(forward)).read()
	lineCountRev = os.popen("cat {} | wc -l".format(reverse)).read()
else:
	logger.warn("{} and {} must both be gzipped or both gunzipped.".format( getIsolateStr(forward, 0), getIsolateStr(reverse, 0)))
	sys.exit(1)	

if(lineCountFwd != lineCountRev):
		logger.warn("Forward and reverse files do not have equal number of reads!")
		sys.exit(1)



#print(lineCountFwd + " " + lineCountRev)

lineIDFwd = re.split(pattern=' ', string=firstLineFwd)
lineIDRev = re.split(pattern=' ', string=firstLineRev)

## print("{} should equal {}".format(lineIDFwd[0], lineIDRev[0]))

fwdR1 = lineIDFwd[1].strip()
revR2 = lineIDRev[1].strip()

##print("R1 has {} and R2 has {}".format(fwdR1, revR2))


## Validate input file run matching: FILTER 3

if(lineIDFwd[0] == lineIDRev[0]):
	logger.info("Input file IDs validated.")
else:
	logger.warn("{} and {} are not paired ends of same sequencing run.".format( getIsolateStr(forward, 0), getIsolateStr(reverse, 0)))
	print("Exiting.")
	sys.exit(1)

## Inspect input file read-pairing: FILTER 4

if(fwdR1.startswith("1") and revR2.startswith("2")):
	logger.info("Input file forward and reverse indicators validated.")
else:
	logger.warn("End-pairing not apparent in {} and {}.".format( getIsolateStr(forward, 0), getIsolateStr(reverse, 0)))



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

##if(re.search(r'\.gz$', forward, flags=re.IGNORECASE)):


outputSuccess = newOutputFolder + "/" + outputFileString
outputFAIL = newOutputFolder + "/" + outputFileString + "FAIL"
runLog = newOutputFolder + "/prinseq.runtime.log"

intMinScore = int(args.min_score)

### Begin performing read clean

logger.info("starting prinseq-lite.pl")

if(re.search(r'\.fastq$', forward, flags=re.IGNORECASE) and re.search(r'\.fastq$', reverse, flags=re.IGNORECASE)):
	if((args.rm_ambig == 'Y') and (args.trim_qual == 'Y')):	
		os.system("prinseq-lite.pl -fastq {} -fastq2 {} -out_good {} -out_bad {} -min_len {} -trim_ns_left 1 -trim_ns_right 1 -ns_max_n {} -trim_qual_right {} -trim_qual_type mean -trim_qual_window 10 -trim_qual_step 5 2> {}"
.format(gunzipForward, gunzipReverse, outputSuccess, outputFAIL,intMinLen, intTrimN, intMinScore, runLog))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'prinseq-lite.pl: reads with > {} Ns removed/trimmed; reads 3-prime end Q < {} trimmed; reads less than {} bp length removed' > {}".format(intTrimN, intMinScore, intMinLen, outputLog))
	elif((args.rm_ambig == 'N') and (args.trim_qual == 'Y')):	
		os.system("prinseq-lite.pl -fastq {} -fastq2 {} -out_good {} -out_bad {} -min_len {} -trim_qual_right {} -trim_qual_type mean -trim_qual_window 10 -trim_qual_step 5 2> {}"
.format(gunzipForward, gunzipReverse, outputSuccess, outputFAIL,intMinLen, intMinScore, runLog))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'prinseq-lite.pl: reads 3-prime end Q < {} trimmed; reads less than {} bp length removed' > {}".format(intMinScore, intTrimN, intMinLen, outputLog))
	elif((args.rm_ambig == 'Y') and (args.trim_qual == 'N')):	
		os.system("prinseq-lite.pl -fastq {} -fastq2 {} -out_good {} -out_bad {} -min_len {} -trim_ns_left 1 -trim_ns_right 1 -ns_max_n {} 2> {}"
.format(forward, reverse, outputSuccess, outputFAIL,intMinLen, intTrimN, runLog))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'prinseq-lite.pl: reads with > {} Ns removed/trimmed; reads less than {} bp length removed' > {}".format(intTrimN, intMinLen, outputLog))
	else:
		os.system("prinseq-lite.pl -fastq {} -fastq2 {} -out_good {} -out_bad {} -min_len {} 2> {}".format(forward, reverse, outputSuccess, outputFAIL, intMinLen, runLog))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'prinseq-lite.pl: reads less than {} bp length removed' > {}".format(intMinLen, outputLog))

elif(re.search('\.fastq.gz', forward, flags=re.IGNORECASE) and re.search('\.fastq.gz', reverse, flags=re.IGNORECASE)):
	gunzipForward = re.sub(r'\.gz$', '', forward)
	gunzipReverse = re.sub(r'\.gz$', '', reverse)

	### Unzip R1
	logger.info("gunzip {}".format(forward))
	os.system("gunzip -c {} > {}".format(forward, gunzipForward))
	logger.info("forward (R1) gunzip complete")
	outputSuccess = newOutputFolder + "/" + outputFileString
	outputFAIL = newOutputFolder + "/" + outputFileString + "FAIL"
	
	### Unzip R2
	logger.info("gunzip {}".format(reverse))
	os.system("gunzip -c {} > {}".format(reverse, gunzipReverse))
	logger.info("reverse (R2) gunzip complete")
	logger.info("forward (R1) and reverse (R2) files trimming start")
	
	### Begin performing read clean
	if((args.rm_ambig == 'Y') and (args.trim_qual == 'Y')):
		os.system("prinseq-lite.pl -fastq {} -fastq2 {} -out_good {} -out_bad {} -min_len {} -trim_ns_left 1 -trim_ns_right 1 -ns_max_n {} -trim_qual_right {} -trim_qual_type mean -trim_qual_window 10 -trim_qual_step 5 2> {}"
.format(gunzipForward, gunzipReverse, outputSuccess, outputFAIL,intMinLen, intTrimN, intMinScore, runLog))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'prinseq-lite.pl: reads with > {} Ns removed/trimmed; reads 3-prime end Q < {} trimmed; reads less than {} bp length removed' > {}".format(intTrimN, intMinScore, intMinLen, outputLog))
	elif((args.rm_ambig == 'N') and (args.trim_qual == 'Y')):
		os.system("prinseq-lite.pl -fastq {} -fastq2 {} -out_good {} -out_bad {} -min_len {} -trim_qual_right {}-trim_qual_type mean -trim_qual_window 10 -trim_qual_step 5 2> {}"
.format(gunzipForward, gunzipReverse, outputSuccess, outputFAIL,intMinLen, intMinScore, runLog))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'prinseq-lite.pl: reads 3-prime end Q < {} trimmed; reads less than {} bp length removed' > {}".format(intMinScore, intTrimN, intMinLen, outputLog))
	elif((args.rm_ambig == 'Y') and (args.trim_qual == 'N')):	
		os.system("prinseq-lite.pl -fastq {} -fastq2 {} -out_good {} -out_bad {} -min_len {} -trim_ns_left 1 -trim_ns_right 1 -ns_max_n {} 2> {}"
.format(gunzipForward, gunzipReverse, outputSuccess, outputFAIL,intMinLen, intTrimN, runLog))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'prinseq-lite.pl: reads with > {} Ns removed/trimmed; reads less than {} bp length removed' > {}".format(intTrimN, intMinLen, outputLog))
	else:
		os.system("prinseq-lite.pl -fastq {} -fastq2 {} -out_good {} -out_bad {} -min_len {} 2> {}".format(gunzipForward, gunzipReverse, outputSuccess, outputFAIL, intMinLen, runLog))
		outputLog = newOutputFolder + "/parameters.log"
		os.system("echo 'prinseq-lite.pl: reads less than {} bp length removed' > {}".format(intMinLen, outputLog))

	os.system("rm -v {}".format(gunzipForward))
	os.system("rm -v {}".format(gunzipReverse))
else:
	logger.warn("Filetype - R1 {} R2 {} must both end in .fastq".format(forward, reverse))
	logger.warn("Filetype - R1 {} and R2 {} must both end in .fastq.gz".format(forward, reverse))
	logger.exception("Input files must be FASTQ ending in either .fastq or .fastq.gz")



if(args.clean_output == 'Y'):
	os.system("rm -v {}'_1_singletons.fastq'".format(outputSuccess))
	os.system("rm -v {}'_2_singletons.fastq'".format(outputSuccess))
	os.system("rm -v {}'FAIL_1.fastq'".format(outputSuccess))
	os.system("rm -v {}'FAIL_2.fastq'".format(outputSuccess))
	os.system("gzip -c {}'_1.fastq' > {}'_R1_001.cleaned.fastq.gz'".format(outputSuccess, outputSuccess))
	os.system("gzip -c {}'_2.fastq' > {}'_R2_001.cleaned.fastq.gz'".format(outputSuccess, outputSuccess))
	os.system("rm -v {}'_1.fastq'".format(outputSuccess))
	os.system("rm -v {}'_2.fastq'".format(outputSuccess))
		
logger.info("prinseq-lite.pl trimming complete")

