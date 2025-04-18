import os, sys, re, csv, statistics
import argparse, logging, warnings
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

## Script for calculating sequence length and average PHRED score per read
## Outputs comma-delimited table to stdout
## Requires Biopython

## Function: A closure for file extension checking

def ext_check(expected_ext, openner):
        def extension(filename):
                if not filename.lower().endswith(expected_ext):
                        raise ValueError()
                return openner(filename)
        return extension

parser = argparse.ArgumentParser(description='Computes sequence lengths and average PHRED for shuffled paired reads in fastq', usage="readShuffledFastq.py filepath/filename.fastq")

parser.add_argument('filename', nargs='+', type=ext_check('.fastq', argparse.FileType('r')))

args = parser.parse_args()

print(args.filename[0].name)

csvRow = []
forwardLen = []
reverseLen = []
forwardAvg = []
reverseAvg = []

iter = 0

myFastq = open(args.filename[0].name, "r")

for record in SeqIO.parse(myFastq, "fastq"):
        #print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
        if(iter % 2 == 0):
                print("%s,%i,%0.2f," % (record.description, len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])), end="")
        elif(iter % 2 == 1):
                if(" " in record.description):
                        strand = record.description.split(" ")
                        print("%s,%i,%0.2f" % (strand[1], len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
                else:
                        print("%i,%0.2f" % ( len(record.seq), statistics.mean(record.letter_annotations["phred_quality"])))
        iter = iter + 1




