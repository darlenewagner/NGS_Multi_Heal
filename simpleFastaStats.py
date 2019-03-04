#!/usr/bin/python

import sys
import os
import os.path
import argparse
import re
import string
import logging
import warnings

## simpleFastaStats.py takes a single fasta-formatted DNA/RNA text file and
## outputs contig count, average contig length, N50 contig lengths, maximum contig length, and cumulative contig length

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

parser = argparse.ArgumentParser(description='output contig count, average contig length, N50 contig, and maximum contig length for input.assembly.fasta', usage="simpleFastaStats.py filepath/input.assembly.fasta --minLength 500(default) --format [brief(default)|verbose|tsv|csv]")

parser.add_argument("filename",type=ext_check('.fasta', argparse.FileType('r')))

parser.add_argument("--minLength", '-min', default='500', type=int)

parser.add_argument("--format", default='brief', type = lambda s : s.lower(), choices=['tsv', 'csv', 'brief', 'verbose', 'verbose_tsv', 'c', 's', 'b', 'v', 'vt'])

args = parser.parse_args()

## call function to extract filename/isolate from inFilePath
inFileName = getIsolateID(args.filename.name)

## Open input file
filehandle = open(args.filename.name, 'r')

intMinLen = args.minLength

draftContigs = []
draftGenome = {}
contigLengths = {}
idCount = 0
contigID = ""
contigStr = "" 

contigCount = 0
contigMax = 0
contigN50 = 0
draftLength = 0
avgContig = 0

for line in filehandle:
	if(re.search(r'^>', line)):
		if(idCount > 0):
			draftGenome[contigID] = contigStr
			contigID = ""
			contigStr = ""
		contigID = line.strip()
		if(re.search(r'\(paired\)', contigID)):
			contigID = contigID.replace('_(paired)', '')
		if(re.search('R1_001_', contigID)):
			contigID = contigID.replace('R1_001_', '')
		draftContigs.append(contigID)
		idCount = idCount + 1
		#print(contigID)
	elif(re.search(r'^(A|T|G|C|U|N)+', line)):
		contigStr = contigStr + line.strip()

draftGenome[contigID] = contigStr

## Close input file
args.filename.close()

for contigKey in draftGenome:
	if( len(draftGenome[contigKey]) > (intMinLen - 1) ):
		contigLengths[contigKey] = len(draftGenome[contigKey])
		##print(contigKey + " => " + str(contigLengths[contigKey]))
	
## new dictionary, sorted from longest contigs to shortest

#contigSortGenome = {}
#sortContigLengths = {}
#sortContigLengths = sorted(contigLengths, reverse=True)

#contigSortGenome = sorted(draftGenome, key=contigLengths.__getitem__, reverse=True)

count = 0

### Loop to find longest contig and contig count given length > intMinLen

for contigID in sorted(contigLengths, key=contigLengths.__getitem__, reverse=True):
	if( contigLengths[contigID] > (intMinLen - 1) ):
		if(count == 0):
			contigMax = contigLengths[contigID]
			top = 1
		count = count + 1
		draftLength = draftLength + contigLengths[contigID]

contigCount = count
avgContig = draftLength/contigCount

### to compute N50, find the contig that 'resides' at 1/2 of draftLength

#contigSortGenome = sorted(draftGenome, key=contigLengths.__getitem__)

cumulativeLength = 0;

for contigID in sorted(contigLengths, key=contigLengths.__getitem__, reverse=True):
	if( contigLengths[contigID] > (intMinLen - 1) ):
		cumulativeLength = cumulativeLength + contigLengths[contigID]
	if(cumulativeLength > (draftLength/2)):
		contigN50 = contigLengths[contigID]
		break

if ( args.format == 'verbose' or args.format == 'v' ):
	print("Assembly File\tMinimum Contig Length:\tcontigCount\tavgContig\tN50\tmaxContig\tdraftLength")
	print("{}\t".format(inFileName), ">", intMinLen - 1 ,"bp:\t", contigCount, "\t", "%.0f" % avgContig, "\t", contigN50, "\t", contigMax, "\t", draftLength)
elif( args.format == 'brief' or args.format == 'b' ):
	print("Assembly\tcontigCount\tavgContig\tN50\tmaxContig")
	print(inFileName + "\t" + str(contigCount) + "\t" + str("%.0f" % avgContig) + "\t" + str(contigN50) + "\t" + str(contigMax))
elif ( args.format == 'tsv' or args.format == 't'):
	print(str(contigCount) + "\t" + str("%.0f" % avgContig) + "\t" + str(contigN50) + "\t" + str(contigMax))
elif ( args.format == 'verbose_tsv' or args.format == 'vt'):
	print("Assembly_File\tcontigCount\tavgContig\tN50\tmaxContig\tdraftLength")
	print("{}\t".format(inFileName), "\t", contigCount, "\t", "%.0f" % avgContig, "\t", contigN50, "\t", contigMax, "\t", draftLength)
elif ( args.format == 'csv' or args.format == 'c' ):
	print(inFileName + "," + str(contigCount) + "," + str("%.0f" % avgContig) + "," + str(contigN50) + "," + str(contigMax))
	
	

