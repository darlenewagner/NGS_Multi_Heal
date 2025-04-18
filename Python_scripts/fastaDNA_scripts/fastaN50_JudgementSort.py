#!/usr/bin/python

import sys
import os
import os.path
import argparse
import re
import string
import logging
import warnings

## fastaJudgementSort.py takes a two fasta-formatted DNA/RNA text files, computes contig counts,
## average contig lengths, N50 contig lengths, and maximum contig lengths for the two files
## the file with the fewest contigs (less fragmented) is retained while the more fragmented
## fasta file is 'damned' to the user-defined --outFile.  If contig count is equal, the file
## with the longer maximium contig is retained

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


parser = argparse.ArgumentParser(description='output contig count, average contig length, N50 contig, and maximum contig length for input.assembly.fasta', usage="fastaJudgementSort.py filepath1/file1.assembly.fasta filepath2/file2.assembly.fasta --outDir hell_path/fasta_hell --minLength 500(default)")

## two input files required
parser.add_argument("filename1",type=ext_check('.fasta', argparse.FileType('r')))

parser.add_argument("filename2",type=ext_check('.fasta', argparse.FileType('r')))


## output folder
parser.add_argument('--outDir', '-D', type=readable_dir, required=True, action='store')

parser.add_argument("--minLength", '-min', default='500', type=int)

parser.add_argument("--format", default='brief', type = lambda s : s.lower(), choices=['mute', 'brief', 'verbose', 'm', 'b', 'v'])

args = parser.parse_args()

## call function to extract filenames/isolates from inFilePath
inFileName1 = getIsolateID(args.filename1.name)
inFileName2 = getIsolateID(args.filename2.name)

helHeim = args.outDir

intMinLen = args.minLength

## Activate logging

logger = logging.getLogger("simpPrinseqLite_R1andR2.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

if(args.filename1.name == args.filename2.name):
	logger.warn("Input files have the same name and path -- trivial comparison!")
	sys.exit(1)

#draftContigs = []

draftGenome1 = {}
contigLengths1 = {}

draftGenome2 = {}
contigLengths2 = {}

idCount1 = 0
idCount2 = 0
contigID = ""
contigStr = "" 

contigCount1 = 0
contigMax1 = 0
contigSecond1 = 0
contigThird1 = 0
contigFourth1 = 0
contigN50_1 = 0
draftLength1 = 0
avgContig1 = 0

contigCount2 = 0
contigMax2 = 0
contigSecond2 = 0
contigThird2 = 0
contigFourth2 = 0
contigN50_2 = 0
draftLength2 = 0
avgContig2 = 0

## Open input file1
filehandle1 = open(args.filename1.name, 'r')

for line in filehandle1:
	if(re.search(r'^>', line)):
		if(idCount1 > 0):
			draftGenome1[contigID] = contigStr
			contigID = ""
			contigStr = ""
		contigID = line.strip()
		if(re.search(r'\(paired\)', contigID)):
			contigID = contigID.replace('_(paired)', '')
		if(re.search('R1_001_', contigID)):
			contigID = contigID.replace('R1_001_', '')
		#draftContigs.append(contigID)
		idCount1 = idCount1 + 1
		#print(contigID)
	elif(re.search(r'^(A|T|G|C|U|N)+', line)):
		contigStr = contigStr + line.strip()

draftGenome1[contigID] = contigStr


## reset contigStr
contigStr = ""

## Close input file1
args.filename1.close()

## Open input file2
filehandle2 = open(args.filename2.name, 'r')

for line in filehandle2:
	if(re.search(r'^>', line)):
		if(idCount2 > 0):
			draftGenome2[contigID] = contigStr
			contigID = ""
			contigStr = ""
		contigID = line.strip()
		if(re.search(r'\(paired\)', contigID)):
			contigID = contigID.replace('_(paired)', '')
		if(re.search('R1_001_', contigID)):
			contigID = contigID.replace('R1_001_', '')
		#draftContigs.append(contigID)
		idCount2 = idCount2 + 1
		#print(contigID)
	elif(re.search(r'^(A|T|G|C|U|N)+', line)):
		contigStr = contigStr + line.strip()

draftGenome2[contigID] = contigStr

## Close input file2
args.filename2.close()


## obtain contig lengths for file1
for contigKey in draftGenome1:
	if( len(draftGenome1[contigKey]) > (intMinLen - 1) ):
		contigLengths1[contigKey] = len(draftGenome1[contigKey])
		##print(contigKey + " => " + str(contigLengths1[contigKey]))

## obtain contig lengths for file2
for contigKey in draftGenome2:
	if( len(draftGenome2[contigKey]) > (intMinLen - 1) ):
		contigLengths2[contigKey] = len(draftGenome2[contigKey])
		##print(contigKey + " => " + str(contigLengths2[contigKey]))
	

##### Obtain contigMax, contigCount, and draftLength

count = 0

  # Loop to find longest file1 contig and contig count given length > intMinLen

for contigID in sorted(contigLengths1, key=contigLengths1.__getitem__, reverse=True):
	if( contigLengths1[contigID] > (intMinLen - 1) ):
		if(count == 0):
			contigMax1 = contigLengths1[contigID]
			top = 1
		if(count == 1):
			contigSecond1 = contigLengths1[contigID]
		if(count == 2):
			contigThird1 = contigLengths1[contigID]
		if(count == 3):
			contigFourth1 = contigLengths1[contigID]
		count = count + 1
		draftLength1 = draftLength1 + contigLengths1[contigID]

contigCount1 = count
avgContig1 = draftLength1/contigCount1

count = 0

  # Loop to find longest file1 contig and contig count given length > intMinLen

for contigID in sorted(contigLengths2, key=contigLengths2.__getitem__, reverse=True):
	if( contigLengths2[contigID] > (intMinLen - 1) ):
		if(count == 0):
			contigMax2 = contigLengths2[contigID]
			top = 1
		if(count == 1):
			contigSecond2 = contigLengths2[contigID]
		if(count == 2):
			contigThird2 = contigLengths2[contigID]
		if(count == 3):
			contigFourth2 = contigLengths2[contigID]
		count = count + 1
		draftLength2 = draftLength2 + contigLengths2[contigID]

contigCount2 = count
avgContig2 = draftLength2/contigCount2


##### Obtain N50

  # to compute N50, find the contig that 'resides' at 1/2 of draftLength

cumulativeLength = 0;

for contigID in sorted(contigLengths1, key=contigLengths1.__getitem__, reverse=True):
	if( contigLengths1[contigID] > (intMinLen - 1) ):
		cumulativeLength = cumulativeLength + contigLengths1[contigID]
	if(cumulativeLength > (draftLength1/2)):
		contigN50_1 = contigLengths1[contigID]
		break

cumulativeLength = 0;

for contigID in sorted(contigLengths2, key=contigLengths2.__getitem__, reverse=True):
	if( contigLengths2[contigID] > (intMinLen - 1) ):
		cumulativeLength = cumulativeLength + contigLengths2[contigID]
	if(cumulativeLength > (draftLength2/2)):
		contigN50_2 = contigLengths2[contigID]
		break


## Output Results unless --format = mute

if ( args.format == 'verbose' or args.format == 'v' ):
	print("Assembly File\tMinimum Contig Length:\tcontigCount\tavgContig\tN50\ttopContig\tsecContig\ttertContig\tdraftLength")
	print("{}\t".format(inFileName1), ">", intMinLen - 1 ,"bp:\t", contigCount1, "\t", "%.0f" % avgContig1, "\t", contigN50_1, "\t", contigMax1, "\t", contigSecond1, "\t", contigThird1, "\t", draftLength1)
	print("{}\t".format(inFileName2), ">", intMinLen - 1 ,"bp:\t", contigCount2, "\t", "%.0f" % avgContig2, "\t", contigN50_2, "\t", contigMax2, "\t", contigSecond2, "\t", contigThird2, "\t", draftLength2)
elif( args.format == 'brief' or args.format == 'b' ):
	print("Assembly\tcontigCount\tavgContig\tN50\tmaxContig")
	print(inFileName1 + "\t" + str(contigCount1) + "\t" + str("%.0f" % avgContig1) + "\t" + str(contigN50_1) + "\t" + str(contigMax1))
	print(inFileName2 + "\t" + str(contigCount2) + "\t" + str("%.0f" % avgContig2) + "\t" + str(contigN50_2) + "\t" + str(contigMax2))

	
### Judgement 

goingToHel = 0

if(contigN50_2 > contigN50_1):
	goingToHel = 1
elif(contigN50_2 < contigN50_1):
	goingToHel = 2
elif(avgContig2 > avgContig1):
	goingToHel = 1
elif(avgContig2 < avgContig1):
	goingToHel = 2
else:
	if(contigMax1 < contigMax2):
		goingToHel = 1
	elif(contigMax1 > contigMax2):
		goingToHel = 2
	elif(contigSecond1 < contigSecond2):
		goingToHel = 1
	elif(contigSecond1 > contigSecond2):
		goingToHel = 2
	elif(contigThird1 < contigThird2):
		goingToHel = 1
	elif(contigThird2 > contigThird2):
		goingToHel = 2
	elif(contigFourth1 < contigFourth2):
		goingToHel = 1
	elif(contigFourth1 > contigFourth2):
		goingToHel = 2

if(goingToHel == 1):
	os.system("mv -v {} {}".format(args.filename1.name, helHeim))
elif(goingToHel == 2):
	os.system("mv -v {} {}".format(args.filename2.name, helHeim))
else:
	os.system("mv -v {} {}".format(args.filename1.name, helHeim))
	
