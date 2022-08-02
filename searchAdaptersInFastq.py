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

## string variable for when an adapter file is not provided
fakeAdapter = "GCCCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

parser = argparse.ArgumentParser(description='Computes sequence lengths and average PHRED for shuffled paired reads in fastq', usage="searchAdaptersInFastq.py filepath/filename.fastq filepath/adapterfile.csv")

parser.add_argument('filename', nargs='+', type=ext_check('.fastq', argparse.FileType('r')))

parser.add_argument('adapterfile', type=ext_check('.csv', argparse.FileType('r')))

args = parser.parse_args()

#print(args.filename[0].name)
adapters = {}
header = []

table = open(args.adapterfile.name, 'r')
adaptable = csv.reader(table, delimiter=',')

matrixAdapters = [row for row in adaptable]

print("Searching for the following adapters:")

idx = 1

while idx < len(matrixAdapters):
        print(matrixAdapters[idx][0], " -> ", matrixAdapters[idx][1])
        idx = idx + 1


myTitle = re.split(r'[\.\/]', args.filename[0].name)

myFastq = open(args.filename[0].name, "r")

idx = 1
iter = 0

searchString = r"\.+" + re.escape(fakeAdapter) + r"$"



for record in SeqIO.parse(myFastq, "fastq"):
        idx = 1
        while idx < len(matrixAdapters):
                if(iter % 2 == 0):
                        if(re.search(r'' + str(matrixAdapters[idx][1]) + '', str(record.seq))):
                                print("%s R1, %s" % (record.id, str(record.seq)))
                elif(iter % 2 == 1):
                        if(re.search(r'' + str(matrixAdapters[idx][1]) + '', str(record.seq))):
                                print("%s R2, %s" % (record.id, str(record.seq)))
                idx = idx + 1
        iter = iter + 1
        
