#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 100GB
#SBATCH --job-name featurecounts
#SBATCH --output /replace/by/path/to/your/scratch/space/SLX-14572/featurecounts/featurecounts.%j.out


featureCounts \
    --primary \
    -C \
    -t exon \
    -g gene_id \
    -a /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.84.gtf \
    -o /replace/by/path/to/your/scratch/space/SLX-14572/featurecounts/SLX-14572.FourSamples.featureCounts \
    /replace/by/path/to/your/scratch/space/SLX-14572/alignment/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam \
    /replace/by/path/to/your/scratch/space/SLX-14572/alignment/SLX-14572.i703_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam \
    /replace/by/path/to/your/scratch/space/SLX-14572/alignment/SLX-14572.i706_i502.HHMJ3BGX3.s_1.r_1.small.fq.bam \
    /replace/by/path/to/your/scratch/space/SLX-14572/alignment/SLX-14572.i706_i503.HHMJ3BGX3.s_1.r_1.small.fq.bam
