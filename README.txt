NGS_Multi_Heal simplifies quality analysis, trimming, and healing of paired-end
reads comprised of forward (R1) and reverse (R2) reads in gzipped fastq files

Prerequisites:
Java 1.7.xxxx or higher
Bpipe 0.9.8.7 or higher
Perl 5.X.X or higher
FASTX-Toolkit-0.0.13 or higher
prinseq 1.20.X or higher

Bpipe download and documentation: http://docs.bpipe.org/
FASTX-Toolkit download and documentation: http://hannonlab.cshl.edu/fastx_toolkit/

### How to run FastX_Trimmer_R1_R2.bpipe with command-line parameters for trim lengths: TRIMR1 and TRIMR2
bpipe run -p TRIMR1=10 -p TRIMR2=20 NGS_Multi_Heal/FastX_Trimmer_R1_R2.bpipe NGS_Multi_Heal/ExampleData/Campy_D5480_R1_001.fastq.gz 
NGS_Multi_Heal/ExampleData/Campy_D5480_R2_001.fastq.gz
## Result: 10 3' bases trimmed from Campy_D5480_R1_001.fastq.gz and 20 3' bases trimmed from Campy_D5480_R2_001.fastq.gz
