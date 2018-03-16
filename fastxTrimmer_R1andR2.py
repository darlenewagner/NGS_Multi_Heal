#!/usr/bin/python

import sys
import os
import argparse
import re

## A command-line script for trimming read pairs by different offsets at 3' ends
## input: filePath/forwardReads.fastq filePath/reverseReads.fastq --trimF int --trimR int
## output: trimmed reads and log files in output folder 'TrimByPython/forwardReads.cleaned.fastq/'
### requires >= FASTX-Toolkit-0.0.13, >= prinseq-lite.pl 0.20.x, and >= perl/5.22.x

outFolder = 'TrimByPython/'

parser = argparse.ArgumentParser(description="trimming forward reads by -t1 and reverse reads by -t2", usage="--trimF / --trimR: require integers from 0 to 20.")

parser.add_argument('--trimF', '-t1', required=True, type=int, choices=range(0,21))
parser.add_argument('--trimR', '-t2', required=True, type=int, choices=range(0,21))

parser.add_argument('forward', type=argparse.FileType('r'))
parser.add_argument('reverse', type=argparse.FileType('r'))

args = parser.parse_args()

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
                isolateString = re.sub(r'R1_001', 'cleaned_fastq', isolateString[0])
        return isolateString

## Simple output file renaming

outputFolder = getIsolateStr(forward)

newOutputFolder = outFolder + outputFolder 

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


if(re.search('\.fastq', forward, flags=re.IGNORECASE) and re.search('\.fastq', reverse, flags=re.IGNORECASE)):
	#print("Nice files for", newOutputFolder)
	
	os.system("fastx_trimmer -Q33")


