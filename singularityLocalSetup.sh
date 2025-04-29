#!/usr/bin/bash

mkdir $PWD/TrimByPython

singularity pull https://depot.galaxyproject.org/singularity/perl:5.32

singularity build my_perl.sif perl\:5.32

singularity pull https://depot.galaxyproject.org/singularity/fastx_toolkit:0.0.14--h503566f_12

singularity build my_fastx_toolkit.sif fastx_toolkit:0.0.14--h503566f_12

singularity pull https://depot.galaxyproject.org/singularity/prinseq:0.20.4--hdfd78af_5

singularity build my_prinseq.sif prinseq:0.20.4--hdfd78af_5
