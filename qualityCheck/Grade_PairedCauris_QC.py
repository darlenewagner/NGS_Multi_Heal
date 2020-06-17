#!/usr/bin/python

import sys
import os.path
import argparse
import re
import csv
from io import StringIO
import pandas as pd
import numpy
#import statistics

## 1. Runs input .fastq.gz file through friends_rachel.pl to find PHRED score and estimated coverage
## 2. Reads fastqc_data.txt of corresponding .fastq.gz file to optain GC%
## 3. print tab-delimited line: file.fastq.gz\tPHRED\tCoverage\tGC%

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

## Function: Check for valid .fastq pair
def properFastq(prospectiveFile):
    if not os.path.isfile(prospectiveFile):
        raise argparse.ArgumentTypeError("prpperFastq:{0} not a valid file".format(prospectiveFile))
    elif not re.match(r'.*_R(1|2)_001\.fastq\.gz', prospectiveFile):
        if not re.match(r'.*_R(1|2)\.fastq\.gz', prospectiveFile):
            raise argparse.ArgumentTypeError("properFastq:{0} must end in _R(1|2)_001.fastq.gz".format(prospectiveFile))
        else:
            return prospectiveFile
    else:
        return prospectiveFile

## Assume R1 and R2 .fastq.gz are in the same folder/filepath

parser = argparse.ArgumentParser(description='Integrate QC results of friends_rachel.pl with fastqc GC output', usage="summary_of_Cauris_QC.py filepath1/WGSisolate_R1_001.fastq.gz filepath1/WGSisolate_R2_001.fastq.gz filepath2/WGSisolate_R1_fastqc/fastqc_data.txt  filepath2/WGSisolate_R2_fastqc/fastqc_data.txt > output.tsv")

## Mutually-exclusive fastQC and --findGC = N
#my_group = parser.add_mutually_exclusive_group(required=True)

## read WGSisolate.fastq.gz file in fasta format
parser.add_argument('forward', type=properFastq, action='store') # type=ext_check('.fastq.gz', argparse.FileType('r')))
parser.add_argument('reverse', type=properFastq, action='store') # type=ext_check('.fastq.gz', argparse.FileType('r')))

## read fastqc_data.txt file
parser.add_argument("fastqcF") #, type=ext_check('.txt', argparse.FileType('r')))
parser.add_argument("fastqcR")

parser.add_argument('--findGC', '-gc', default='Y', choices=['Y', 'N'], help="parse fastQC output files (Y|N)?")

parser.add_argument('--format', '-f', default='tsv', choices=['tsv', 't', 'csv', 'c'], help="Tab-delimited or comma-delimited output (tsv|csv)?")

parser.add_argument('--stdout', '-s', default='N', choices=['Y', 'N'], help="Write to file or STDOUT (Y|N)?")

args = parser.parse_args()

forward = args.forward
reverse = args.reverse

isolateIDF = getIsolateID(forward)
isolateIDR = getIsolateID(reverse)

isolateOut = re.sub("(_L001)?_R1(_001)?", "", isolateIDF)

if((args.format == 'csv') or (args.format == 'c')):
    isolateOut = isolateOut + '.csv'
elif((args.format == 'tsv') or (args.format == 't')):
    isolateOut = isolateOut + '.tsv'

#print(isolateOut)

fastqc1 = args.fastqcF
fastqc2 = args.fastqcR

#if(args.findGC == 'Y'):
gcPercent1 = os.popen("head -10 {} | tail -1".format(fastqc1)).read()
gcPercent2 = os.popen("head -10 {} | tail -1".format(fastqc2)).read()
simpleGC1 = gcPercent1.split()
simpleGC2 = gcPercent2.split()
simpGC1 = simpleGC1[1]
simpGC2 = simpleGC2[1]


readMetricsF = os.popen("perl /scicomp/groups/OID/NCEZID/DFWED/MDB/data_analysis/ydn3/ROSS/scripts/friends_rachel.pl --fast -e 12400000 {} | tail -1".format(forward)).read()
readMetricsR = os.popen("perl /scicomp/groups/OID/NCEZID/DFWED/MDB/data_analysis/ydn3/ROSS/scripts/friends_rachel.pl --fast -e 12400000 {} | tail -1".format(reverse)).read()

simpReadMetricsF = readMetricsF.split()
simpReadMetricsR = readMetricsR.split()

coverageF = simpReadMetricsF[7]
coverageR = simpReadMetricsR[7]

coverF = float(coverageF)
coverR = float(coverageR)

phredF = simpReadMetricsF[5]
phredR = simpReadMetricsR[5]

phredFstatus = 'FAIL'
phredRstatus = 'FAIL'
coverageStatus = 'FAIL'
simpGC1status = 'FAIL'
simpGC2status = 'FAIL'

if(float(phredF) > 28.0):
    phredFstatus = 'PASS'
if(float(phredR) > 20.0):
    phredRstatus = 'PASS'
if((round(coverF, 0) + round(coverR, 0)) > 19.99):
    coverageStatus = 'PASS'
if((42 < int(simpGC1)) and (int(simpGC1) < 47) and (args.findGC == 'Y')):
    simpGC1status = 'PASS'
if((42 < int(simpGC2)) and (int(simpGC2) < 47) and (args.findGC == 'Y')):
    simpGC2status = 'PASS'



with open(isolateOut, 'w') as outfile:
    if(args.stdout == 'N'):
        sys.stdout = outfile
    if((args.format == 'tsv') or (args.format == 't')):
        if(args.findGC == 'Y'):
            print("WGS isolate\tPHRED\tPHRED Status\tCoverage\tCoverage Status\tG+C%\tG+C% Status")
            print(isolateIDF + "\t" + phredF + "\t" + phredFstatus + "\t" + str(round(coverF, 0)) + "x\t" + coverageStatus + "\t" + simpGC1 + "%" + "\t" + simpGC1status)
            print(isolateIDR + "\t" + phredR + "\t" + phredRstatus + "\t" + str(round(coverR, 0)) + "x\t" + coverageStatus + "\t" + simpGC2 + "%" + "\t" + simpGC2status)
        else:
            print("WGS isolate\tPHRED\tPHRED Status\tCoverage\tCoverage Status")
            print(isolateIDF + "\t" + phredF + "\t" + phredFstatus + "\t" + str(round(coverF, 0)) + "x\t" + coverageStatus)
            print(isolateIDR + "\t" + phredR + "\t" + phredRstatus + "\t" + str(round(coverR, 0)) + "x\t" + coverageStatus)
    else:
        if(args.findGC == 'Y'):
            print("WGS isolate, PHRED, PHRED Status, Coverage, Coverage Status, G+C%, G+C% Status")
            print(isolateIDF + ", " + phredF + ", " + phredFstatus + ", " + str(round(coverF, 0)) + "x, " + coverageStatus + ", " + simpGC1 + "%" + ", " + simpGC1status)
            print(isolateIDR + ", " + phredR + ", " + phredRstatus + ", " + str(round(coverR, 0)) + "x, " + coverageStatus + ", " + simpGC2 + "%" + ", " + simpGC2status)
        else:
            print("WGS isolate, PHRED, PHRED Status, Coverage, Coverage Status")
            print(isolateIDF + ", " + phredF + ", " + phredFstatus + ", " + str(round(coverF, 0)) + "x, " + coverageStatus)
            print(isolateIDR + ", " + phredR + ", " + phredRstatus + ", " + str(round(coverR, 0)) + "x, " + coverageStatus)

outfile.close()
