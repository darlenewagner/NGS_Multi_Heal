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

parser.add_argument("--stdout", '-s', default='N', choices=['Y', 'N'] )

args = parser.parse_args()

inFilePath = args.filename.name

## call function to extract filename/isolate from inFilePath
inFileName = getIsolateID(inFilePath)

outFileName = re.sub("(\.sorted)?\.bam", "", inFileName)

if(args.format == 'verbose' or args.format == 'v' or args.format == 'tsv' or args.format == 't'):
        outFileName = outFileName + '.tsv'
elif(args.format == 'csv' or args.format == 'c'):
        outFileName = outFileName + '.csv'
else:
        outFileName = outFileName + '.txt'


sumCommand = "samtools view -h {} | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{{sum+=$1}}END{{print sum}}'"

refLength = os.popen(sumCommand.format(inFilePath)).read()

percMapString = os.popen("samtools flagstat {} | head -5 | tail -1".format(inFilePath)).read()

percMapArr = re.split('(\(|\s:)', percMapString)
percMapFloat = float(percMapArr[2][:5])

#print(percMapArr[2])
#print(percMapFloat)

intRefLength = int(float(refLength))

avgDepth = os.popen("samtools depth {} | awk '{{sum+=$3}}END{{print sum / {} }}'".format(inFilePath, intRefLength)).read()

fAvgDepth = float(avgDepth)

grade = 'FAIL'
if((percMapFloat > 49.99) and (fAvgDepth > 19.50)):
    grade = 'PASS'

with open(outFileName, 'w') as outfile:
    if(args.stdout == 'N'):
        sys.stdout = outfile
    if((args.format == 'verbose') or (args.format == 'v')):
        print("isolate\treference length\treference coverage\tpercent reads mapped\tstatus")
        print(inFileName + "\t" + "%.0ibp" % intRefLength + "\t" + "%.2fx" % fAvgDepth + "\t" + "%.2f" % percMapFloat + "%" + "\t" + grade)
    elif((args.format == 'tsv') or (args.format == 't')):
        print("isolate\treference coverage\tpercent reads mapped\tstatus")
        print(inFileName + "\t" + "%.2fx" % fAvgDepth + "\t" + "%.2f" % percMapFloat + "%" + "\t" + grade)
    elif((args.format == 'csv') or (args.format == 'c')):
        print("isolate, reference coverage, percent reads mapped, status")
        print(inFileName + ", " + "%.2fx" % fAvgDepth + ", " + "%.2f" % percMapFloat + "%, " + grade)
    elif((args.format == 'simp') or (args.format == 's')):
        print("{}\t".format(inFileName), "%.2f" % fAvgDepth)
outfile.close()

args.filename.close()

