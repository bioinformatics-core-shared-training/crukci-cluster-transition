#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 512
#SBATCH --job-name hello_world
#SBATCH --output /replace/by/path/to/your/scratch/space/hello_world.%j.out

echo Hello World from the cluster!
