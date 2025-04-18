#!/usr/bin/python

import sys
import os
import os.path
import argparse
import re
import logging
import warnings

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

parser = argparse.ArgumentParser(description='Find contig count, N50, and max contig for input.assembly.fasta', usage="simpleFastaStats.py filepath/input.assembly.fasta --format [simp(default)|verbose]")

parser.add_argument("filename",type=ext_check('.fasta', '.fa', argparse.FileType('r')))

parser.add_argument("--format", default='simp', type = lambda s : s.lower(), choices=['simp', 'verbose', 's', 'v'])

args = parser.parse_args()

draftGenome = []

with open(filename, 'r') as genome:
	for line in genome:
		if()




