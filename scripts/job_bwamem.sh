#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 100GB
#SBATCH --job-name bwamem
#SBATCH --output /replace/by/path/to/your/scratch/space/SLX-14572/alignment/bwamem.i706_i517.small.%j.out

/home/bioinformatics/software/bwa/bwa-0.7.15/bwa mem -M -t 4 -R "@RG\tID:1\tLB:SLX-14572.i706_i517\tSM:SLX-14572\tPU:HHMJ3BGX3.1" \
    /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38/bwa/hsa.GRCh38 \
    /replace/by/path/to/your/scratch/space/SLX-14572/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.gz \
    | /home/bioinformatics/software/samtools/samtools-1.6/bin/samtools view -b \
    > /replace/by/path/to/your/scratch/space/SLX-14572/alignment/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam
