# Session 2: Cluster

## Learning Objectives

- Get familiar with most useful shell commands from session 1
- View sequencing data and its associated file formats
- Rename files using loop
- Run commands in a shell script
- Submit jobs to the cluster using SLURM
- Run quality control, artefact removal and alignment on the cluster


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

### File formats

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


## Tracking your work - always create a `README.txt` file

- Create a `README.txt` file in each of your working directory
- Describe what you are trying to achieve
- Copy/paste the commands you ran using your `history`

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data.
> - Create a `README.txt` file.
> - Type a small description and the commands you used to retrieve the data.
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Renaming multiple files - always keep your original data

We've seen in session 1 how to use wildcard `*` to get a command executed on multiple files, as well as combining them using pipe `|`. Unfortunately it is not possible to use wildcard to rename files. We will have to use a loop to do some operation once for each thing in a list.

> :computer: **EXERCISE** Go to your Terminal window, or open a new one.
>
> - Let's get back to Nelle's data from session 1. On the command line download the [zipped data file](https://github.com/bioinformatics-core-shared-training/crukci-cluster-transition/raw/master/session1-data.zip) using `wget` followed by `unzip` to decompress the archive and navigate to the `session1-data/nelle/creatures` directory using `cd` command and list all the files using `ls`.
>
> :tada: Congratulations! :thumbsup: You did it! :wink:

There are two files in this directory but imagine we had hundred which we would like to rename to `original_*.bat` to take a back up of them before editing them.

If you run

```
mv *.dat original_*.dat
```
you'll get an message about how to use `mv` because `mv` cannot receives more than two inputs. Instead, we can use a loop to do some operation once for each thing in a list. Let's start by printing the name of each of us using the `echo` command. We could do it one command at a time, or we can use a loop `for variable_name in element_1 element_2 element_3; do command $variable_name; done` to repeat the same command three times for each element of the list:
```
echo Anne
echo Rob
echo Jochen
echo Katie
echo Ummi

for name in Anne Rob Jochen Katie Ummi; do echo $name; done
```

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and go to `session1-data/nelle/creatures`.
>
> - Write a loop to print out the name of each file first and then,
> - Write a loop to rename these files `original_*.bat`
>
> :tada: Congratulations! :thumbsup: You did it! :wink:

It is really important to name your variable with a **meaningful name** and not a random one or a one letter one. Programs are only useful if people can understand them, so meaningless names (like `x`) or misleading names (like `temperature`) increase the odds that the program won't do what its readers think it does.


## Simple shell script

We are finally ready to see what makes the shell such a powerful programming environment. We are going to take the commands we repeat frequently and save them a file so that we can re-run all those operations again later by typing a single command. For historical reasons, a bunch of commands saved in a file is usually called a **shell script**, but make no mistake: these are actually small programs.

Let's go back to the `session1-data/nelle/molecules` to extract lines 11 to 15 of each PDB file using a shell script called `middle.sh`.

First, we have to create the file `middle.sh`
```shell
cd session1-data/nelle/molecules
nano middle.sh
```
and type the commands we want to run
```shell
head -15 octane.pdb | tail -5
```
make this file executable to you by changing its mode
```shell
chmod u+x middle.sh
```
and finally run the command
```shell
./middle.sh
```

What if we want to select lines from an arbitrary file? We could edit `middle.sh` each time to change the filename, but that would probably take longer than just retyping the command. Instead, let's edit `middle.sh` and replace `octane.pdb` with a very special variable called `$1`. `$1` means the first parameter on the command line. We can now run our script like this:
```shell
./middle.sh octane.pdb
```

We still need to edit `middle.sh` each time we want to adjust the range of lines, though. Let's fix that by using the special variables `$2` and `$3`.

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and go to `session1-data/nelle/molecules`.
>
> - Update the script `middle.sh` to take two other parameters on the command line for the range of lines to select
> - Run the script
>
> :tada: Congratulations! :thumbsup: You did it! :wink:

This works, but it may take the next person who reads `middle.sh` a moment to figure out what it does. We can improve our script by adding some **comments** at the top.

```shell
# Select lines from the middle of a file.
# Usage: middle.sh filename -end_line -num_lines
head $2 $1 | tail $3
```


## Cluster job submission using SLURM, a job scheduler

### What is SLURM (Simple Lightweight Unix Resource Manager)?

Slurm is an open source, fault-tolerant, and highly scalable cluster management and job scheduling system for large and small Linux clusters.

We are mostly interested by its job scheduling aspect which allocate resources to a user for a specified amount of time. Slurm provides resource management for the processors allocated to a job, so that multiple job steps can be simultaneously submitted and queued until there are available resources within the job's allocation.

### Key SLURM commands

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

### Submit a job

We are going to start by submitting a very simple job to the cluster using the `echo` command and a shell script `job.sh` containing specific Slurm instructions.

First, log in onto the cluster head node
```shell
ssh my_username@clust1-headnode.cri.camres.org
```
create a `job.sh` file containing these lines, do replace `/scratcha/xxlab/my_username` by your scratch space
```shell
#!/bin/sh
#SBATCH --partition general
#SBATCH --mem 512
#SBATCH --job-name hello_world
#SBATCH --output /scratcha/xxlab/my_username/hello_world.%j.out

echo Running on $(hostname): 'Hello world!'
```
and finally submit your job to the cluster using the command `sbatch`
```shell
sbatch job.sh
```

In the `job.sh` script, we have specific `SBATCH` instructions:
- `--partition`: select a specific queue to submit the job
- `--mem`: specify the real memory required e.g. `2GB` is `2048`
- `--job-name`: specify a name for the job allocation
- `--output`: connect the batch script's standard output directly to the file name specified. By default both standard output and standard error are directed to the same file. The filename is `hello_world.%j.out` where the `%j` is replaced by the job ID.

See [sbatch](http://slurm.schedmd.com/sbatch.html) man page for all the options and explanation on submitting a batch script to Slurm.

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and go to `session1-data/nelle/`.
>
> - Copy `molecules/` onto your scratch space on the cluster using `scp -r`
> - Submit `middle.sh` script to Slurm to extract lines 20-23 of `octane.pdb`
> - Write a loop to submit jobs for all PDB files by modifying `job.sh` to take one command line argument which will be given at each step of the loop
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Quality control

We are now going back to our own data and check read quality using FastQC. Here we will use the command line version of FastQC to check what kind of options this tool has. This command will display all the parameters you can use when running `FastQC`.
```shell
/home/bioinformatics/software/fastqc/fastqc-v0.11.5/fastqc --help
```
Now we will run a very simple command:
```shell
/home/bioinformatics/software/fastqc/fastqc-v0.11.5/fastqc -o /scratcha/xxlab/my_username/SLX-ID/fastqc/ --noextract -f fastq my_file.fastq.gz
```

The options are:
- `-o FOLDER_NAME`: you can define this way in which folder you want FastQC to output its results
- `–noextrac`t: tells FastQC not to uncompress the output file after creating it
- `-f fastq`: defines that the input file is a fastq file (valid formats are bam, sam, bam_mapped, sam_mapped and fastq)

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data.
> - Create a `job.sh` to run FastQC
> - Send job to the cluster and wait for the results
> - Update your `README.txt` file with what you've done
> - View the html report in a web browser, you may have to copy back this file on your own computer to be able to view it using the `scp` command
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Artefact removal

Once you had a closer look at the quality report you may realise that the data quality is not too bad, however we still might be able to improve it with a quality based trimming since the quality usually drops towards the end of the reads. We will use `Cutadapt` for trimming.

As the current Bioinformatics' core installation is broken, the easiest is to install it in your home directory on the cluster:
```shell
pip install --user --upgrade cutadapt
```
Let's have a look at the help page
```shell
~/.local/bin/cutadapt --help
```
In some case, all we want to do is to remove low quality bases from our reads. We can use the following command to do this:
```shell
~/.local/bin/cutadapt -m 10 -q 20 -o my_file_trimmed.fastq.gz my_file.fastq.gz
```

The parameters are:
- `-m 10`: discards all reads that would be shorter than a read length of 10 after the trimming
- `-q 20`: trims low-quality bases from the 3’ end of the reads; if two comma-separated cutoffs are given, the 5’ end is trimmed with the first cutoff, the 3’ end with the second
- `-o FILE_NAME`: the output file name

Once the trimming has finished we will want to check the quality of our trimmed reads as well to make sure, we are happy with its results: the trimming improved the quality and it didn’t introduce new artefacts.

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data.
> - Create a new `job.sh` to run Cutadapt
> - Send job to the cluster and wait for the results
> - Update your `README.txt` file with what you've done
> - Rerun FastQC and view the new html report in a web browser
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Sequencing alignment using BWA

```shell
/home/bioinformatics/software/bwa/bwa-0.7.15/bwa
```

```shell
/home/bioinformatics/software/bwa/bwa-0.7.15/bwa mem -M -t 4 /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38/bwa/hsa.GRCh38.bwt my_file_trimmed.fastq.gz > my_file_trimmed.fastq.sam
```


## Reference materials

- [CRUK Bioinformatics Autumn School 2017: Functional Genomics](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/)
  - [Reference genomes and common file formats, Dóra Bihary, Sept17](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/Introduction/SS_DB/Materials/Lectures/Lecture1_fileFormats_DB.pdf)
  - [Quality control and artefact removal](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/Introduction/SS_DB/Materials/Lectures/Lecture2_qualityControl_artefactRemoval_DB.pdf)
