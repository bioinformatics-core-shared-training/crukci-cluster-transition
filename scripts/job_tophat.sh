#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 100GB
#SBATCH --job-name tophat
#SBATCH --output /replace/by/path/to/your/scratch/space/SLX-14572/alignment_tophat/tophat.i706_i517.small.%j.out

/home/bioinformatics/software/tophat/tophat-2.1.1/tophat \
  --GTF /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/annotation/hsa.hg38.gtf \
  --bowtie1 --min-anchor "3" --num-threads "4" \
  --tmp-dir /replace/by/path/to/your/scratch/space/SLX-14572/alignment_tophat/tmp/ \
  --output /replace/by/path/to/your/scratch/space/SLX-14572/alignment_tophat/ \
  /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/tophat/hsa.hg38  \
  /replace/by/path/to/your/scratch/space/SLX-14572/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.gz
