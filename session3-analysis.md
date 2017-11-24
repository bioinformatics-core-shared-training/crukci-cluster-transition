# Session 3: Analysis steps by steps


## Learning Objectives

- Submit jobs to the cluster using SLURM
- Get cluster file systems mounted on your local macOS
- Run quality control, artefact removal and alignment on the cluster


## Submit a job

We did submit a very simple job to the cluster using the `echo` command and a shell script `job.sh` containing specific Slurm instructions. Please check the instructions on [Can I submit jobs onto the cluster?](can-i-submit-jobs.md) and run it again to make sure you can submit jobs to the cluster and understand how to do it.

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


## Get cluster file systems mounted on your macOS

[FUSE for macOS](https://osxfuse.github.io/) allows you to extend macOS's native file handling capabilities via third-party file systems. Combining with [SSHFS](https://github.com/osxfuse/osxfuse/wiki/SSHFS), it gets file system accessible via `ssh` mounted directly on your macOS computer.

Install:
- Latest release of [FUSE for macOS](https://github.com/osxfuse/osxfuse/releases) by downloading the `.dmg` file
- Latest release of [SSHFS](https://github.com/osxfuse/sshfs/releases) by downloading the `.pkg` file

Open a Terminal window, and create a mount point in your home directory `mnt/scratchb` for example:
```
cd ~                    # go to home directory
mkdir -p mnt/scratchb   # create intermediate directories as required with -p option
cd mnt/scratchb         # go to this directory
pwd                     # return current working directory name
```

Mount your cluster `/mnt/scratchb/my_group/my_username` onto your local machine using:
```
sshfs my_username@clust1-headnode.cri.camres.org:/mnt/scratchb/my_group/my_username /Users/my_username/mnt/scratchb
```
To list all mounted file system, use:
```
mount
```
To unmount `scratchb`, use:
```
umount
```

Aliases can be created in your `~/.profile` to save you time and to avoid typing complex commands every time you need them. You could enter these lines at the beginning of your `~/.profile` file using your preferred editor [atom](https://atom.io/) which could be launch from the command line by typing `atom`:
```
atom ~/.profile
```
Add enter these lines into the file:
```
### Aliases
alias mntclustsb='sshfs pajon01@clust1-headnode.cri.camres.org:/mnt/scratchb/my_group/my_username /Users/my_username/mnt/scratchb'
alias umntclustsb='umount /Users/my_username/mnt/scratchb'
```
Save and open a new Terminal window, your new aliases are now available as 'new' commands! :thumbsup:


## FastQ quality control

We are now going back to our own data and check read quality using FastQC. Here we will use the command line version of FastQC to check what kind of options this tool has. This command will display all the parameters you can use when running `FastQC`.
```
/home/bioinformatics/software/fastqc/fastqc-v0.11.5/fastqc --help
```
Now to run FastQC on your downloaded sequencing data in `SLX-ID` folder, you will need to run this command but we are not running this command directly, we will submit this job to the cluster:
```
/home/bioinformatics/software/fastqc/fastqc-v0.11.5/fastqc -o /scratcha/xxlab/my_username/SLX-ID/ --noextract -f fastq my_file.fastq.gz
```
If you need to get sequencing data again, check the steps from [Session 1: Shell - Getting sequencing: Using CRUKCI infrastructure data](session1-shell.md#using-crukci-infrastructure)

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
> - View the html report in a web browser, you may have to copy back this file on your own computer to be able to view it using the `scp` command or mount your scratch space using `sshfs`
>
> :tada: Congratulations! :thumbsup: You did it! :wink:

When running jobs on the cluster, you may wish to keep an eye on the output by using `tail -f` which appends data as the file grows:
```
tail -f my_output_file_name.out
```
Check your job is still running and in the queue using [squeue](http://slurm.schedmd.com/squeue.html), you can combine it with `grep my_username` to only extract information about your jobs:
```
squeue | grep my_username
```
You can display information of all your submitted jobs using [sacct](http://slurm.schedmd.com/sacct.html) but most importantly you need to check that your job has completed using:
```
sacct -j JobID
```
If you wish to kill your job before it completes, run [scancel](http://slurm.schedmd.com/scancel.html) using:
```
scancel JobID
```


## Artefact removal in reads

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
> - Create a new `job_cutadapt.sh` to run Cutadapt
> - Send job to the cluster
> - Check the output while running
> - Update your `README.txt` file with what you've done
> - Rerun FastQC and view the new html report in a web browser and let's spot the differences
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Short sequencing read alignment to a reference genome

Let's prepare a small fastQ file to avoid waiting for too long for the alignment to run.

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data.
> - Combine `zcat` and `head` to extract 100,000 reads from a fastQ file
> - Redirect the output into a new filename e.g. `SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq`
> - Compress your output using `gzip`
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


### ChipSeq - alignment with BWA

[BWA](http://bio-bwa.sourceforge.net/bwa.shtml) can map sequences against a large reference genome, such as the human genome. BWA MEM can map longer sequences (70bp to 1Mbp) and is generally recommended for high-quality queries as it is faster and more accurate.

Our installed version is located in `/home/bioinformatics/software/bwa/bwa-0.7.15/bwa`, to get specific help from the command line for BWA MEM run `/home/bioinformatics/software/bwa/bwa-0.7.15/bwa mem` or for general help `/home/bioinformatics/software/bwa/bwa-0.7.15/bwa`.

For 75 bp reads against GRCh38 reference genome, we are going to run `bwa mem`:
```
/home/bioinformatics/software/bwa/bwa-0.7.15/bwa mem -M -t 4 -R "@RG\tID:1\tLB:SLX-14572.i706_i517\tSM:SLX-14572\tPU:HHMJ3BGX3.1" \
    /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38/bwa/hsa.GRCh38 \
    /scratchb/xxlab/my_username/SLX-14572/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.gz \
    > /scratchb/xxlab/my_username/SLX-14572/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.sam
```
- `-M`: leaves the best (longest) alignment for a read as is but marks additional alignments for the read as secondary
- `-t 4` number of processor cores
- `-R "@RG\tID:1\tLB:SLX-14572.i706_i517\tSM:SLX-14572\tPU:HHMJ3BGX3.1"` add read group header to identify your aligned reads which will help when merging bam files later, `ID` for identifier, `LB` for library identifier (SLX-ID) and barcode, `SM` for sample name and `PU` for platform unit including flow cell and lane number

The output of BWA is a SAM file, [samtools](http://www.htslib.org/) to convert it to bam format, we will be using `samtools view -b`, our installed version is located here `/home/bioinformatics/software/samtools/samtools-1.6/bin/samtools`. We are going to pipe the output of `bwa mem` into `samtools` to avoid writing multiple files on disk:
```
/home/bioinformatics/software/bwa/bwa-0.7.15/bwa mem -R @RG\tID:1\tLB:SLX-14572.i706_i517\tSM:SLX-14572\tPU:HHMJ3BGX3.1 -t 4 \    
    /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38/bwa/hsa.GRCh38 \
    /scratchb/xxlab/my_username/SLX-14572/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.gz \
    | /home/bioinformatics/software/samtools/samtools-1.6/bin/samtools view -b \
    > /scratchb/xxlab/my_username/SLX-14572/alignment/SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam
```

To view the header of your aligned reads, you can use `samtools view -H my_file.bam` or to view some aligned reads use `samtools view my_file.bam | tail -10`.

If you have short sequence reads (< 70bp), you will need to run two steps `bwa aln` followed by `bwa samse` for single end or `bwa sampe` for paired end data.

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data.
> - Create an `alignment` directory
> - Create a new `job.sh` to run bwa mem on the small fastQ file created earlier `SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq`
> - Send job to the cluster
> - Check the output while running
> - Update your `README.txt` file with what you've done
>
> :tada: Congratulations! :thumbsup: You did it! :wink:

If you're getting failed job, you may wish to check that you have enough disk space on your scratch space using:
```
lfs quota -h /scratchb/xxlab/
```

### RNASeq - alignment with TopHat

- See [TopHat](http://ccb.jhu.edu/software/tophat/index.shtml) documentation.
- Our installed version is located here `/home/bioinformatics/software/tophat/tophat-2.1.1/tophat`
- Running TopHat:
  ```
  tophat \
      --GTF $GeneAnnotationFile \
      --bowtie1 --min-anchor "3" --num-threads "4" \
      --tmp-dir $TempDirectory \
      --output $OutputDirectory \
      $TopHatIndexPrefix  \
      $FastqFile

  $GeneAnnotationFile=
  /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/annotation/hsa.hg38.gtf

  $TopHatIndexPrefix=
  /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/tophat/hsa.hg38
  ```

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data.
> - Create an new `alignment_tophat` directory
> - Create a new `job.sh` to run tophat on the small fastQ file `SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq`
> - Send job to the cluster
> - Check the output while running
> - Update your `README.txt` file with what you've done while job is running
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Viewing alignment data

- Aligned sequence formats: SAM/BAM/CRAM
- SAM format stands for Sequence Alignment Map format and it is the standard format for aligned sequence data which is recognised by the majority of software and browsers

Explore alignment data files from [SAM/BAM and related specifications](http://samtools.github.io/hts-specs/) and the [SAM Format Specification](http://samtools.github.io/hts-specs/SAMv1.pdf).

To view the header of your aligned reads, you can use `samtools view -H my_file.bam` or to view some aligned reads use `samtools view my_file.bam | tail -10`.

- See [samtools](http://www.htslib.org/) document for more information.
- Our installed version is located here `/home/bioinformatics/software/samtools/samtools-1.6/bin/samtools`.

## Alignment quality metrics, sort and mark duplicates using Picard tools

- See [Picard tools](http://broadinstitute.github.io/picard/) documentation.
- Our installed version is located here `/home/bioinformatics/software/picard/picard-2.14.0/picard.jar`, you need `java` installed
- Running Picard tools:
  - Alignment Metrics
  ```
  java -jar PICARDJAR CollectAlignmentSummaryMetrics \
      I=InputBam \
      O=OutputMetricsFile \
      REF_FLAT=GenomeReferenceInFASTA
  ```
  - InsertSize Metrics
  ```
  java -jar PICARDJAR CollectInsertSizeMetrics \
      I=InputBam \
      O=OutputMetricsFile \
      H=OutputHistogramPlotPDF \
      VALIDATION_STRINGENCY=SILENT
  ```
  - Sort bam
  ```
  java -jar PICARDJAR SortSam \
      I=InputBam \
      O=OutputBam \
      SORT_ORDER=coordinate - One of {unsorted, queryname, coordinate, duplicate, unknown}
  ```
  - Mark Duplicates (keep or delete)
  ```
  java -jar PICARDJAR MarkDuplicates \
      I=InputBam \
      O=OutputBam \
      M=OutputMetricsFile \
      REMOVE_DUPLICATES=false
  ```
  - RNAseq Metrics
  ```
  java -jar PICARDJAR CollectRnaSeqMetrics \
      I=InputBam \
      O=OutputMetricsFile \
      REF_FLAT=GenomeReferenceInFASTA \
      STRAND_SPECIFICITY="NONE" \
      ASSUME_SORTED="false"  \ (true is quicker)
      VALIDATION_STRINGENCY=SILENT
  ```
  - WGS metrics
  ```
  java -jar PICARDJAR CollectWgsMetrics \
      I=InputBam \
      O=OutputMetricsFile \
      REF_FLAT=GenomeReferenceInFASTA
  ```

## Take home message: everyday cluster commands

```
sbatch job.sh
sacct -j JobID
tail -f my_output_file_name.JobID.out
```

## Reference materials

- [CRUK Bioinformatics Autumn School 2017: Functional Genomics](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/)
  - [Quality control and artefact removal](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/Introduction/SS_DB/Materials/Lectures/Lecture2_qualityControl_artefactRemoval_DB.pdf)
  - [Short Read Alignment](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/Introduction/SS_DB/Materials/Lectures/Lecture3_ShortRead_Alignment_SS.pdf)
