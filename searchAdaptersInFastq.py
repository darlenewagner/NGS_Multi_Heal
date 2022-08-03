import os, sys, re, csv, statistics
import argparse, logging, warnings
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
import pandas as pd
#from collections import defaultdict

def ext_check(expected_ext, openner):
        def extension(filename):
                if not filename.lower().endswith(expected_ext):
                        raise ValueError()
                return openner(filename)
        return extension

parser = argparse.ArgumentParser(description='Searches for adapters in fastq shuffled paired reads. Requires Bio.SeqIO.', usage="searchAdaptersInFastq.py filepath/filename.fastq filepath/adapterfile.csv")

parser.add_argument('filename', type=ext_check('.fastq', argparse.FileType('r')))

parser.add_argument('adapterfile', type=ext_check('.csv', argparse.FileType('r')))

parser.add_argument('--outputType', '-o', default='C', choices=['C', 'R'], help="--outputType C for counts of adapter occurrences and R for reads containing adapters.")

parser.add_argument('--verbose', '-v', default='N', choices=['Y', 'N'], help="--verbose enables output of adapters list at beginning.")

args = parser.parse_args()

myTitle = re.split(r'\/', args.filename.name)
newTitle = re.sub('\.f(ast)?q(\.gz)', '', myTitle[len(myTitle) - 1])

adapters = {}
header = []

table = open(args.adapterfile.name, 'r')
adaptable = csv.reader(table, delimiter=',')

matrixAdapters = [row for row in adaptable]

if(args.verbose == 'Y'):
        print("Searching for the following adapters:")

idx = 1

countAdapters = {}

while idx < len(matrixAdapters):
        if(args.verbose == 'Y'):
                print(matrixAdapters[idx][0], " -> ", matrixAdapters[idx][1])
        countAdapters[str(matrixAdapters[idx][0])] = 0
        idx = idx + 1




myFastq = open(args.filename.name, "r")

idx = 1
iter = 0

# to quantify forward (R1), reverse (R2) strand bias
countStrand = {'R1' : 0, 'R2' : 0}

readsWithAdapters = open('/scicomp/home-pure/ydn3/test_Python3.9.1/test_Biopython/' + newTitle + '.with_adapters.fastq', 'w')

for record in SeqIO.parse(myFastq, "fastq"):
        idx = 1
        while idx < len(matrixAdapters):
                strand = record.description.split(" ")
                if(re.search(r'^1', strand[1])):
                        if(re.search(r'' + str(matrixAdapters[idx][1]) + '', str(record.seq))):
                                if(args.outputType == 'R'):
                                        ## Save record on read with adapter to reads.with_adapters.fastq
                                        countStrand['R1'] = countStrand['R1'] + 1
                                        # print(record.format("fastq"), end="")
                                        readsWithAdapters.write(record.format("fastq"))
                                elif(args.outputType == 'C'):
                                        countAdapters[str(matrixAdapters[idx][0])] = countAdapters[str(matrixAdapters[idx][0])] + 1
                elif(re.search(r'^2', strand[1])):
                        if(re.search(r'' + str(matrixAdapters[idx][1]) + '', str(record.seq))):
                                if(args.outputType == 'R'):
                                        countStrand['R2'] = countStrand['R2'] + 1
                                        # print(record.format("fastq"), end="")
                                        readsWithAdapters.write(record.format("fastq"))
                                elif(args.outputType == 'C'):
                                        countAdapters[str(matrixAdapters[idx][0])] = countAdapters[str(matrixAdapters[idx][0])] + 1
                idx = idx + 1
        iter = iter + 1
        
readsWithAdapters.close()

if(args.outputType == 'R'):
        ## output strand bias to standard output
        print("Forward-Reverse Read Bias:")
        print("%s: %i, %s: %i" % ("R1", countStrand['R1'], "R2", countStrand['R2']))
if(args.outputType == 'C'):
        ## output adapter occurrence counts to standard output
        print("Adapter_Name, Count")
        for k in countAdapters.keys():
                print(k + ", " + str(countAdapters[k]))
