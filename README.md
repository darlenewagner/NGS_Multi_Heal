## NGS_Multi_Heal
NGS_Multi_Heal simplifies quality analysis, trimming, and healing of paired-end reads
comprised of forward (R1) and reverse (R2) reads in gzipped or unzipped fastq files.

## Preprocessing with NGS_Multi_Heal
* Trimming 3-prime ends of reads prior to high-quality Single Nucleotide Polymorphisms (hqSNPs) analysis
* Removing ambiguous nucleotides (Ns) from reads prior to hqSNPs analysis or de novo assembly
* Removing reads under a user-specified minimum prior to de novo assembly

## Prerequisites
* Perl 5.X.X or higher
* Python 3.X.X or higher
* Java 1.7.xxxx or higher
* Bpipe 0.9.8.7:  [bpipe.org](http://docs.bpipe.org)
* FASTX-Toolkit-0.0.13: [hannonlab.cshl.edu](http://hannonlab.cshl.edu/fastx_toolkit)
* prinseq 1.20.X: [sourceforge.net/projects/prinseq](https://sourceforge.net/projects/prinseq/files/standalone/)
* CG-Pipeline Perl: [github.com/lskatz](https://github.com/lskatz/CG-Pipeline)

## Example Preprocessing for hqSNPs
##### Trim 3' ends of forward and reverse reads by 5 and 15 base pair positions, respectively:
```python fastxTrimmer_R1andR2.py reads_R1_001.fastq reads_R2_001.fastq --trimF 5 --trimR 15 -outDir trimmed/```
##### Remove reads < 40 bp and containing > 0 ambiguous nucleotides:
```python simpPrinseqLite_R1andR2.py reads_R1_001.fastq reads_R2_001.fastq --min_len 40 --rm_ambig Y --ambig_allow 0 --outDir trimmed/```
## Example Preprocessing for de novo assembly
##### Remove reads < 100 bp, containing > 0 ambiguous nucleotides, and trim regions with quality < 24:
```python qualPrinseqLite_R1andR2.py reads_R1_001.fastq reads_R2_001.fastq --min_len 100 --rm_ambig Y --ambig_allow 0 --trim_qual Y --min_score 24 --outDir trimmed/```

## Coming Soon
* Python wrapper for SPAdes BayesHammer
* Python scripts for managing DNA .fasta or .txt
