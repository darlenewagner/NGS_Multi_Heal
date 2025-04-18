## simpleFastaStats.py takes a single fasta-formatted DNA/RNA text file and
## outputs contig count, average contig length, N50 contig lengths, maximum contig length, and cumulative contig length

## fastaJudgementSort.py takes a two fasta-formatted DNA/RNA text files, computes contig counts,
## average contig lengths, N50 contig lengths, and maximum contig lengths for the two files
## the file with the fewest contigs (less fragmented) is retained while the more fragmented
## fasta file is 'damned' to the user-defined --outFile.  If contig count is equal, the file
## with the longer maximium contig is retained

