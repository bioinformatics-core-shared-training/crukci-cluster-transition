#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 512
#SBATCH --job-name fastqc
#SBATCH --output /replace/by/path/to/your/scratch/space/SLX-14572/fastqc.i706_i517.%j.out

/home/bioinformatics/software/fastqc/fastqc-v0.11.5/fastqc \
  -o /replace/by/path/to/your/scratch/space/SLX-14572/fastqc/ \
  --noextract \
  -f fastq \
  /replace/by/path/to/your/scratch/space/SLX-14572/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.fq.gz
