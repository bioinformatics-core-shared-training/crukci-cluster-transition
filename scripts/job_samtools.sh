#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 100GB
#SBATCH --job-name samtools
#SBATCH --output /replace/by/path/to/your/scratch/space/SLX-14572/alignment/samtools.i706_i517.small.%j.out

/home/bioinformatics/software/samtools/samtools-1.6/bin/samtools sort \
  -T sorted \
  /replace/by/path/to/your/scratch/space/SLX-14572/alignment/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam \
  > /replace/by/path/to/your/scratch/space/SLX-14572/alignment/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq_sorted.bam; \
  /home/bioinformatics/software/samtools/samtools-1.6/bin/samtools index \
  /replace/by/path/to/your/scratch/space/SLX-14572/alignment/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq_sorted.bam
