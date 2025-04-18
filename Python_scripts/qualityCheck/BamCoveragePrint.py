#!/usr/bin/python

import sys
import os.path
import argparse
import re

## Function: A closure for file extension checking

def ext_check(expected_ext, openner):
	def extension(filename):
		if not filename.lower().endswith(expected_ext):
			raise ValueError()
		return openner(filename)
	return extension

## Function: Filename extractor from filepath
def getIsolateID(filePathString):
	splitStr = re.split(pattern='/', string=filePathString)
	fileNameIdx = len(splitStr) - 1
	isolateString = re.split(pattern='\.', string=splitStr[fileNameIdx]) 
	if(len(isolateString[0]) < 10):
		isolateString = re.split(pattern='\.', string=splitStr[0])
	return isolateString[0]


parser = argparse.ArgumentParser(description='Show coverage of reads mapped in BAM file', usage="BamCoveragePrint.py filepath/file.bam --format [simp(default)|csv|tsv|verbose]")

## Specify positional parameter, 'filename', and invoke function ext_check to validate .bam extension
parser.add_argument("filename",type=ext_check('.bam', argparse.FileType('r')))

parser.add_argument("--format", default='simp', type = lambda s : s.lower(), choices=['simp', 'csv', 'tsv', 'verbose', 's', 'c', 't', 'v'] )

args = parser.parse_args()

inFilePath = args.filename.name

## call function to extract filename/isolate from inFilePath
inFileName = getIsolateID(inFilePath)

sumCommand = "samtools view -h {} | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{{sum+=$1}}END{{print sum}}'"

refLength = os.popen(sumCommand.format(inFilePath)).read()

intRefLength = int(float(refLength))

avgDepth = os.popen("samtools depth {} | awk '{{sum+=$3}}END{{print sum / {} }}'".format(inFilePath, intRefLength)).read()

fAvgDepth = float(avgDepth)

# verbose output
if ( args.format == 'verbose' or args.format == 'v' ):
	print(inFileName, "\tBAM Length: ", "%.0ibp" % intRefLength, "\tBAM Coverage: ", "%.2fx" % fAvgDepth)

# output as .tsv
elif ( args.format == 'tsv' or args.format == 't'):
	print(inFileName, "\t", "%.0ibp" % intRefLength, "\t", "%.2fx" % fAvgDepth)

# output as .csv
elif ( args.format == 'csv' or args.format == 'c' ):
	print("{},".format(inFileName), "{}bp,".format(intRefLength), "%.2fx" % fAvgDepth)

# output coverage only
elif ( args.format == 'simp' or args.format == 's' ):
	print("{}\t".format(inFileName), "%.2f" % fAvgDepth)


args.filename.close()

