#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 512
#SBATCH --job-name cutadapt
#SBATCH --output /scratcha/xxlab/my_username/SLX-14572/cutadapt.i706_i517.%j.out

~/.local/bin/cutadapt -m 10 -q 20 -o my_file_trimmed.fastq.gz my_file.fastq.gz
