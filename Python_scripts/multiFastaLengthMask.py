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

parser = argparse.ArgumentParser(description='output contig count, average contig length, N50 contig, and maximum contig length for multiple input.assembly.fasta', usage="multiFastaStats.py filepath/input.assembly.fasta --minLength 500(default) --format [brief(default)|verbose|tsv|csv]")

parser.add_argument("filename",type=ext_check('.fasta', argparse.FileType('r')))

parser.add_argument("--minLength", '-min', default='200', type=int)

parser.add_argument("--format", default='brief', type = lambda s : s.lower(), choices=['tsv', 'csv', 'brief', 'verbose', 'c', 's', 'b', 'v'])

## arrays of dict type variables
#GenomeDrafts = []
#GenomeContigs = []

inFileName = []

draftContigs = []
draftGenome = {}
contigLengths = {}
idCount = 0
contigID = ""
contigStr = "" 

contigCount = []
maxContig = []
contigN50 = []
drLength = 0
draftLength = []
avgContig = []

args = parser.parse_args()

intMinLen = args.minLength

basePath = os.getcwd()
print(args.filename)
inPath = os.path.split(args.filename.name)
filehandle = open(args.filename.name, 'r')
baseName = os.path.basename(args.filename.name)
baseBaseName = os.path.splitext(baseName)[0]

idxFile = 0

##### Begin first loop, read input file #####

for line in filehandle:
        if(re.search(r'^>', line)):
                if(idCount > 0):
                        draftGenome[contigID] = contigStr
                        contigID = ""
                        contigStr = ""

                contigID = line.strip()
                draftContigs.append(contigID)
                idCount = idCount + 1
        elif(re.search(r'^(A|T|G|C|U|N)+', line)):
                contigStr = contigStr + line.strip()


        

draftGenome[contigID] = contigStr

## End first inner loop
## Close input file

filehandle.close()	

### Second loop to populate dict of contig lengths and output contigs over --minLength

for contigKey in draftGenome:
        if( len(draftGenome[contigKey]) > (intMinLen - 1) ):
                contigLengths[contigKey] = len(draftGenome[contigKey])
                tempFile = re.sub(">", "", contigKey)
                tempFileName = tempFile.split('_cov')
                print("{}".format(tempFileName[0]))
                outFasta = open(basePath + '/' + inPath[0] + '/' + tempFileName[0] + '.fasta', 'w')
                outFasta.write( ">{}_{}\n".format(tempFileName[0], baseBaseName) )
                outFasta.write( "{}\n".format(draftGenome[contigKey]) )
                outFasta.close()

### End second loop


