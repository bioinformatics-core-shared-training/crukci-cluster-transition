# Session 2: Cluster

## Learning Objectives

- Get familiar with most useful shell commands from session 1
- Submit job to the cluster using SLURM
-

## File formats

- Reference genomes: consist of a mixture of known chromosomes and unplaced contigs called Genome Reference Assembly
  - Path to reference genomes: `/scratchb/bioinformatics/reference_data/reference_genomes/`
  - Path to assembly: `/scratchb/bioinformatics/reference_data/reference_genomes/$organism/$assembly` e.g. for Human GRCh38 `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38`

- Unaligned sequences: Fasta and FastQ
  - FastQ: sequence information with [quality scores](https://en.wikipedia.org/wiki/FASTQ_format)

- Aligned sequences: SAM/BAM/CRAM
- Summarised genomic features
  - Genomic intervals: BED
  - Gene annotation: GFF/GTF
  - Genomic scores: Wiggle files, BEDgraphs, BigWigs

## Viewing sequencing data

```shell
ssh my_username@clust1-headnode.cri.camres.org           # access the cluster head node
cd /scratchb/my_group/my_username/                       # go to your scratch space
java -jar /home/my_username/clarity-tools.jar -l SLX-ID  # download your sequencing data
cd SLX-ID/fastq/                                         # navigate to fastq data folder
zcat my_sequence_file.fq.gz | more                       # output the content of the file paging through text one screenful at a time
```

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your scratch folder.
> - Download your preferred project data.
> - Visualise one FASTQ file
>
> :tada: Congratulations! :thumbsup: You did it! :wink:

FASTQ file consists of multiple blocks of these four lines
```
@NS500222:320:HHMJ3BGX3:1:11101:24390:1371 1:N:0:TAAGGCGA+ATAGAGAG
CATCTGCAAGTTGGAGACCCAGATAAGCCAGTAATGTAGTTCAGTCCATGACCAAACTGTCTCTTATACACATCT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEE
```
- Line 1: `@` followed by identifier
- Line 2: Sequence information
- Line 3: Character `+`
- Line 4: [Quality scores](https://en.wikipedia.org/wiki/FASTQ_format)

## Tracking your work - best practices

- Create a `README.txt` file in each of your working directory
- Explain what you are trying to do
- Copy the most important commands you've ran using your `history`

## Renaming multiple files - always keep your original data

We've seen in session 1 how to use wildcard `*` to get a command executed on multiple files, as well as combining them using pipe `|`. Unfortunately it is not possible to use wildcard to rename files. We will have to use a loop to do some operation once for each thing in a list.

> :computer: **EXERCISE** Go to your Terminal window, or open a new one.
>
> - Let's get back Nelle's data from session 1, on the command line download the [zipped data file](https://github.com/bioinformatics-core-shared-training/crukci-cluster-transition/raw/master/session1-data.zip) using `wget` followed by `unzip` to decompress the archive and navigate to the `session1-data/nelle/creatures` directory using `cd` command and list all the files using `ls`.
>
> :tada: Congratulations! :thumbsup: You did it! :wink:

There are two files in this directory but imagine we had hundred which we would like to rename to `original_*.bat` to take a back up of them before editing them.

If you run

```
mv *.dat original_*.dat
```
you'll get an message about how to use `mv` because `mv` cannot receives more than two inputs. Instead, we can use a loop to do some operation once for each thing in a list. Let's start by printing the name of each file using the `echo` command of each file in turn:
```
for filename in basilisk.dat unicorn.dat; do echo $filename; done
```

## Simple shell script


## Cluster job submission using SLURM, a job scheduler

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

Create a `job.sh` file
```shell
#!/bin/sh
#SBATCH --no-requeue
#SBATCH -p general
#SBATCH -J 000000000-BC2Y2_alignment_pipeline
#SBATCH --mem 4096
#SBATCH --mincpus 1
#SBATCH --open-mode truncate
#SBATCH -o /mnt/scratcha/bioinformatics/solexa/Runs/171110_M01712_0413_000000000-BC2Y2/alignment/pipeline.%j.out

# autoanalysis generated shell script
export MEM_VALUE=4096
export MEM_LIMIT=$[${MEM_VALUE}*1024]
export JAVA_OPTS="-Xmx$[${MEM_VALUE}-128]M -Xms$[${MEM_VALUE}-128]M"

/home/solexa/sequencingpipelines//alignment/bin/run-pipeline --mode=slurm /mnt/scratcha/bioinformatics/solexa/Runs/171110_M01712_0413_000000000-BC2Y2/alignment/run-meta.xml
```

Submit job to the cluster
```shell
sbatch job.sh
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



## Quality control

## Artefact removal

## Sequencing alignment with BWA

## Reference materials

- [CRUK Bioinformatics Autumn School 2017: Functional Genomics](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/)
  - [Reference genomes and common file formats, Dóra Bihary, Sept17](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/Introduction/SS_DB/Materials/Lectures/Lecture1_fileFormats_DB.pdf)
  - [Quality control and artefact removal](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/Introduction/SS_DB/Materials/Lectures/Lecture2_qualityControl_artefactRemoval_DB.pdf)
