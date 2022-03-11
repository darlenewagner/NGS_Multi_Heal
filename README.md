
## NGS_Multi_Heal
NGS_Multi_Heal simplifies quality analysis, trimming, and healing of paired-end reads
comprised of forward (R1) and reverse (R2) reads in gzipped or unzipped fastq files.  
Two scripts in this package, fastxTrimmer_R1andR2.py and windowQualPrinseqLite_R1andR2.py,
were employed in read-healing analysis in Wagner et al. (2021) PeerJ. 9:e12446 
(https://doi.org/10.7717/peerj.12446)

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

## Example Preprocessing of reads for hqSNPs
### Method I. Uniform 3-prime trim for forward and reverse reads 
##### Trim 3' ends of forward and reverse reads by 5 and variable (XX) base pair positions, respectively:
```python fastxTrimmer_R1andR2.py reads_R1_001.fastq reads_R2_001.fastq --trimF 5 --trimR XX -outDir trimmed/```
* trim forward (R1) reads by 5 bp.
* trim reverse (R2) reads by 5 bp when average PHRED quality score is > 31.00
* otherwise, trim reverse reads by 10 bp for PHRED qual. < 31.00
* trim reverse reads by 15 bp for PHRED qual. < 30.00
* trim reverse reads by 20 bp for PHRED qual. < 29.00
### Method II. 
##### Step 1. Remove reads < 40 bp and/or containing > 0 ambiguous nucleotides:
```python simpPrinseqLite_R1andR2.py reads_R1_001.fastq reads_R2_001.fastq --min_len 40 --rm_ambig Y --ambig_allow 0 --outDir trimmed/```
##### Step 2. Trim 3' ends of forward and reverse reads as shown under Method I:
```python fastxTrimmer_R1andR2.py trimmed/reads_prinseq_R1_001.fastq trimmed/reads_prinseq_R2_001.fastq --trimF 5 --trimR XX -outDir trimmed/```
## Example Preprocessing genomeA reads for de novo assembly
### Method I. Single master wrapper for prinseq-lite.pl 
##### Remove reads < 100 bp, containing > 0 ambiguous nucleotides, and trim regions with quality < 26:
```python qualPrinseqLite_R1andR2.py genomeA_reads_R1_001.fastq genomeA_reads_R2_001.fastq --min_len 100 --rm_ambig Y --ambig_allow 0 --trim_qual Y --min_score 26 --outDir trim_genomeA_reads/```
### Method II. Prinseq-lite.pl wrapper followed by fastx_trimmer wrapper
##### Step 1. Remove reads < 100 bp and/or containing > 1 ambiguous nucleotides
```python simpPrinseqLite_R1andR2.py genomeA_reads_R1_001.fastq genomeA_reads_R2_001.fastq --min_len 100 --rm_ambig Y --ambig_allow 1 --outDir trim_genomeA_reads/```
##### Step 2. Trim 3' ends of forward and reverse reads as shown under Method I for hqSNPs
```python fastxQualAdaptTrimmer_R1andR2.py trim_genomeA_reads/genomeA_reads_R1_001_prinseq/genomeA_reads_R1_001_prinseq_1.fastq trim_genomeA_reads/genomeA_reads_R1_001_prinseq/genomeA_reads_R1_001_prinseq_2.fastq --trimF 5 --trimR XX --outDir trim_genomeA_reads/```
## Coming Soon
* Python wrapper for SPAdes BayesHammer
* Python scripts for managing DNA .fasta or .txt
