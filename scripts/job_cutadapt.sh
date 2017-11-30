#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 512
#SBATCH --job-name cutadapt
#SBATCH --output /replace/by/path/to/your/scratch/space/SLX-14572/cutadapt.i706_i517.%j.out

~/.local/bin/cutadapt -m 10 -q 20 -o my_file_trimmed.fastq.gz my_file.fastq.gz
