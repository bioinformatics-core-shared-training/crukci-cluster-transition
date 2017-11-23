#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 512
#SBATCH --job-name middle_pdb
#SBATCH --output /scratcha/xxlab/my_username/molecules/middle_pdb.%j.out

./middle.sh $1 23 4
