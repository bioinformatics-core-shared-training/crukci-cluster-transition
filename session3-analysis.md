# Session 3: Analysis steps by steps


## Learning Objectives

- Submit jobs to the cluster using SLURM
- Run quality control, artefact removal and alignment on the cluster


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

- [CRUK Bioinformatics Autumn School 2017: Functional Genomics](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/)
  - [Quality control and artefact removal](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/Introduction/SS_DB/Materials/Lectures/Lecture2_qualityControl_artefactRemoval_DB.pdf)
