#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 512
#SBATCH --job-name fastqc
#SBATCH --output /scratcha/xxlab/my_username/SLX-14572/fastqc.i706_i517.%j.out

/home/bioinformatics/software/fastqc/fastqc-v0.11.5/fastqc -o /scratcha/xxlab/my_username/SLX-14572 --noextract -f fastq SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.fq.gz
