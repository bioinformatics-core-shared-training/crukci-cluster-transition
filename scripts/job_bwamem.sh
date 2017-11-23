#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 100GB
#SBATCH --job-name bwamem
#SBATCH --output /scratchb/xxlab/my_username/SLX-14572/alignment/bwamem.i706_i517.small.%j.out

/home/bioinformatics/software/bwa/bwa-0.7.15/bwa mem -M -t 4 \
  /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38/bwa/hsa.GRCh38 \
  /scratchb/xxlab/my_username/SLX-14572/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.gz \
  > /scratchb/xxlab/my_username/SLX-14572/alignment/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.sam
