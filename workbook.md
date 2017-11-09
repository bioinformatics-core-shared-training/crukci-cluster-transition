# How-to use the CRUK-CI cluster

Transitioning from Galaxy to Cluster usage.
This page hold information about how to use CRUK-CI 2016 cluster.

## Training material on the intranet of CRUK-CI

http://bioinf-wiki001/doku.php?id=newstarters:computationalinfrastructure:clusterone

## Connecting to Cluster

You first need a cluster account, request one via Helpdesk - IT <ithelpdesk@cruk.cam.ac.uk>.

```shell
ssh clust1-headnode.cri.camres.org
```

## Cluster directories' structure

- Home directory: `/home/$username`
- Shared installed software: `/home/bioinformatics/software`
- Scratch/working areas: `/mnt/scratcha` and `/mnt/scratchb`
- Reference data: `/mnt/scratchb/bioinformatics/reference_data`

The two scratch/working areas are separate but equivalent. These areas are (deliberately) not backed, they are a massively parallel distributed file system
called Lustre. Pick one to do your work, and to make best use vary which one to use.

Large files may need to be stripped to improve performance or files that many jobs reads. See `lfs setstripe --help` for help.

Limit number of files in directory, if possible avoid 10,000s files in single directory.

## Reference genomes

- Path to reference genomes: `/scratchb/bioinformatics/reference_data/reference_genomes/`
- Path to assembly: `/scratchb/bioinformatics/reference_data/reference_genomes/$organism/$assembly` e.g. for Human GRCh38 `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38`

What Bioinformatics Core maintains:
- Genome sequence (fasta)
- Alignment indices: BWA, TopHat, Bowtie (1, 2)
- Annotations:
  - GTF format gene model
  - RefFlat format gene model
  - Signal artifact list (if available)

## Getting sequencing data

The Bioinformatics Core provides a tool for downloading files for projects, libraries and runs that you can use from the command line or integrate into your Java application. This is available internally from:
http://intranet.cri.camres.org/core-facilities/bioinformatics/sequencing/api

Save [this file](http://internal-bioinformatics.cruk.cam.ac.uk/software/clarity-tools.jar) to your working area. You can run the tool from the command line with:

```shell
wget http://internal-bioinformatics.cruk.cam.ac.uk/software/clarity-tools.jar
java -jar clarity-tools.jar --help
```

Usage example: `java -jar clarity-tools.jar -l SLX-14572`

Another tool developed by the Bioinformatics Core is [kickstart](http://intranet.cri.camres.org/core-facilities/bioinformatics/sequencing/kickstart) which is  designed for downstream analysis instead of just downloading.

```shell
/home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart --help
```

Usage example: `/home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart -l SLX-14572`

## Using the scheduler for job submission

### The most useful SLURM (Simple Lightweight Unix Resource Manager) commands

Key commands for the impatient:
- [sbatch](http://slurm.schedmd.com/sbatch.html): Submit jobs
- [squeue](http://slurm.schedmd.com/squeue.html): View the queue
- [sacct](http://slurm.schedmd.com/sacct.html): View jobs' state information
- [scancel](http://slurm.schedmd.com/scancel.html): Kill jobs


Type command followed by -h for command line usage details.

```shell
sbatch -h
squeue -h
sacct -h
scancel -h
```

[SLURM command summary](https://slurm.schedmd.com/pdfs/summary.pdf)

### Simple parallel

Solving many similar and independent tasks:
- Analysis split into tasks
- Task assigned to one cpu
- No inter-task communication
- More throughput by running more tasks
- Task runtime varies

90% of bioinformatics codes fall into this model

### Submit a job

TODO

## Typical work flow

```shell
ssh clust1-headnode.cri.camres.org
cd /scratchb/bioinformatics/pajon01

# Getting sequencing data
mkdir SLX-14572
cd SLX-14572
/home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart -l SLX-14572

# To get data for running alignment pipeline for example using BWA-MEM
/home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart -l SLX-14572 -a bwamem -s homo_sapiens -v grch38
```
