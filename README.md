## NGS_Multi_Heal
NGS_Multi_Heal simplifies quality analysis, trimming, and healing of paired-end reads
comprised of forward (R1) and reverse (R2) reads in gzipped or unzipped fastq files.

## Prerequisites
* Perl 5.X.X or higher
* Python 3.X.X or higher
* Java 1.7.xxxx or higher
* Bpipe 0.9.8.7:  [bpipe.org](http://docs.bpipe.org)
* FASTX-Toolkit-0.0.13: [hannonlab.cshl.edu](http://hannonlab.cshl.edu/fastx_toolkit)
* prinseq 1.20.X: [sourceforge.net/projects/prinseq](https://sourceforge.net/projects/prinseq/files/standalone/)
* CG-Pipeline Perl: [github.com/lskatz](https://github.com/lskatz/CG-Pipeline)

## Examples
##### Trim 3' ends of forward and reverse reads by XX and YY base pair positions, respectively:
```python fastxTrimmer_R1andR2.py reads_R1_001.fastq reads_R2_001.fastq --trimF XX --trimR XX -outDir trimmedReads/```
##### Parse .bam file for coverage of reads to reference - output as .csv to terminal:
```python BamCoveragePrint.py reads_to_ref.sorted.bam --format csv```
##### Find percentage of reads with ambiguous nucleotides (N):
```python countReadsWithAmbig.py NGSreads_R1_001.fastq --format percent```

## Coming Soon
* Python wrapper for prinseq-lite
* Python wrapper for SPAdes BayesHammer
* Python scripts for managing DNA .fasta or .txt
