#!/usr/bin/python

import sys
import os
import argparse
import re

## A command-line script for trimming read pairs by different offsets at 3' ends
## input: filePath/forwardReads.fastq filePath/reverseReads.fastq --trimF int --trimR int
## output: trimmed reads and log files in output folder 'forwardReads.cleaned.fastq/'
### requires >= FASTX-Toolkit-0.0.13, >= prinseq-lite.pl 0.20.x, >= perl/5.22.x, and >= Bpipe 0.9.8.x

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

#intTrimFwd = int(float(trimFwd))
#intTrimRev = int(float(trimRev))

# try ARGPARSE, so that trim args are flags, not positional parameters, a numcpus flag can also be used
# ARGPARSE also has a --help Usage statement and --outDir to allow user to name their own output folder

os.system("bpipe run -p TRIMR1={} -p TRIMR2={} /scicomp/home/ydn3/NGS_Multi_Heal/FastX_Trimmer_R1_R2.bpipe {} {}".format(intTrimFwd, intTrimRev, forward, reverse))

# change filename, forward, to folder name used by bpipe run

pathName = forward.split('/')

lenName = len(pathName)

folderName = re.sub(pattern = 'fastq.gz', string = pathName[lenName - 1], repl = 'cleaned.fastq')

#print(folderName)

os.system("bpipe log > /scicomp/home/ydn3/NGS_Multi_Heal/{}/fastxTrimBpipe.log".format(folderName))


