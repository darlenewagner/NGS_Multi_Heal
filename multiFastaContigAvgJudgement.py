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

## Function: Checks existence of --outDir
def readable_dir(prospective_dir):
	if not os.path.isdir(prospective_dir):
    		raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
	if os.access(prospective_dir, os.R_OK):
		if( not prospective_dir.endswith("/") ):
			prospective_dir = prospective_dir + "/"
		return prospective_dir
	else:
		raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))


parser = argparse.ArgumentParser(description='compare average contig and contig counts among multiple .fasta, move lower quality assemblies to Hel', usage="multiFastaContigAvgJudgement.py filepath/input.assembly*.fasta --minLength 500(default) --format [brief(default)|verbose|tsv|csv]")

parser.add_argument("filename",type=ext_check('.fasta', argparse.FileType('r')), nargs='+')

## output folder
parser.add_argument('--outDir', '-D', type=readable_dir, required=True, action='store')

## minimum contig length
parser.add_argument("--minLength", '-min', default='500', type=int)

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

helHeim = args.outDir

intMinLen = args.minLength

idxFile = 0

##### Begin logging #####

logger = logging.getLogger("multiFastaContigAvgJudgement.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)

if(len(args.filename) < 2):
	print("Input Error: Two or more .fasta files required!")
	sys.exit(1)

##### Begin multiple input file loop #####

for filehandle in args.filename:
	inFileName.append(getIsolateID(filehandle.name))
	
	## First inner loop to read input file lines

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

	### End first inner loop

	## Close input file
	filehandle.close()

	### Second inner loop to populate dict of contig lengths

	for contigKey in draftGenome:
		if( len(draftGenome[contigKey]) > (intMinLen - 1) ):
			contigLengths[contigKey] = len(draftGenome[contigKey])
			##print(contigKey + " => " + str(contigLengths[contigKey]))

	### End second innner loop

	count = 0

	### Third inner loop to find longest contig and contig count given length > intMinLen

	for contigID in sorted(contigLengths, key=contigLengths.__getitem__, reverse=True):
		if( contigLengths[contigID] > (intMinLen - 1) ):
			if(count == 0):
				maxContig.append(contigLengths[contigID])
				top = 1
			count = count + 1
			drLength = drLength + contigLengths[contigID]
	draftLength.append(drLength)
	
	### End third inner loop
       
	contigCount.append(count)

	avgContig.append(draftLength[idxFile]/contigCount[idxFile])

	### to compute N50, find the contig that 'resides' at 1/2 of draftLength
	
	drLength = 0
	cumulativeLength = 0;

	### Fourth inner loop to calculate N50
	
	for contigID in sorted(contigLengths, key=contigLengths.__getitem__, reverse=True):
		if( contigLengths[contigID] > (intMinLen - 1) ):
			cumulativeLength = cumulativeLength + contigLengths[contigID]
		if(cumulativeLength > (draftLength[idxFile]/2)):
			contigN50.append(contigLengths[contigID])
			break
	
	### End fourth inner loop

	draftContigs = []
	draftGenome = {}
	contigLengths = {}
	idCount = 0
	contigID = ""
	contigStr = "" 

	idxFile = idxFile + 1

##### End of multiple input file loop #####	

idx = 0

for idx in range(len(inFileName)):
	if ( args.format == 'verbose' or args.format == 'v' ):
		print("Assembly File\tMinimum Contig Length:\tcontigCount\tavgContig\tN50\tmaxContig\tdraftLength")
		print("{}\t".format(inFileName[idx]), ">", intMinLen - 1 ,"bp:\t", contigCount[idx], "\t", "%.0f" % avgContig[idx], "\t", contigN50[idx], "\t", maxContig[idx], "\t", draftLength[idx])
	elif( args.format == 'brief' or args.format == 'b' ):
		print("Assembly\tcontigCount\tavgContig\tN50\tmaxContig")
		print(inFileName[idx] + "\t" + str(contigCount[idx]) + "\t" + str("%.0f" % avgContig[idx]) + "\t" + str(contigN50[idx]) + "\t" + str(maxContig[idx]))
	elif ( args.format == 'tsv' or args.format == 't'):
		print(str(contigCount[idx]) + "\t" + str("%.0f" % avgContig[idx]) + "\t" + str(contigN50[idx]) + "\t" + str(maxContig[idx]))
	elif ( args.format == 'csv' or args.format == 'c' ):
		print(inFileName[idx] + "," + str(contigCount[idx]) + "," + str("%.0f" % avgContig[idx]) + "," + str(contigN50[idx]) + "," + str(maxContig[idx]))
	idx = idx + 1

	
### Judgement Day ### 
### All assembly files are sent to Hel except for theWorthy ###

## index of the file with optimal assembly metrics
theWorthy = 0

for helCount in range( 1, len(inFileName) ):
	if(avgContig[helCount] > avgContig[theWorthy]):
		theWorthy = helCount
		##print(theWorthy, " ", avgContig[helCount], " ", avgContig[theWorthy])
	elif(contigCount[helCount] < contigCount[theWorthy]):
		theWorthy = helCount
			
idx = 0

for idx in range( len(inFileName) ):
	if(idx != theWorthy):
		os.system("mv -v {} {}".format(args.filename[idx].name, helHeim))

