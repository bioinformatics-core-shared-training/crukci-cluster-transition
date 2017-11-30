#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 512
#SBATCH --job-name middle_pdb
#SBATCH --output /replace/by/path/to/your/scratch/space/molecules/middle_pdb.%j.out

./middle.sh $1 23 4
