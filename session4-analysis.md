# Session 4: Analysis steps by steps

## Learning Objectives

- Run analysis steps on the cluster including
  - Artefact removal in reads
  - Alignment
  - Sort
  - Mark duplicates
  - Alignment quality metrics
- View alignment data in IGV

## Artefact removal in reads

Once you had a closer look at the quality report you may realise that the data quality is not too bad, however we still might be able to improve it with a quality based trimming since the quality usually drops towards the end of the reads. We will use `Cutadapt` for trimming.

As the current Bioinformatics' core installation is broken, the easiest is to install it in your home directory on the cluster:
```
pip install --user --upgrade cutadapt
```
Let's have a look at the help page
```
~/.local/bin/cutadapt --help
```
In some case, all we want to do is to remove low quality bases from our reads. We can use the following command to do this:
```
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
> - Check the output while running and your job using `sacct` on the cluster
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

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node. Keep two terminal windows open onto the cluster when possible.
>
> - Navigate to your project data.
> - Create an `alignment` directory
> - Create a new `job.sh` to run bwa mem on the small fastQ file created earlier `SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq`
> - Send job to the cluster
> - Check the output while running and your job using `sacct` on the cluster
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
- Before running it you will need to add `bowtie` onto your path, by adding three symbolic links from `bowtie` installation directory to your `~/bin/` using the command `ln -s` (It does work because `~/bin/` is on our path, try `echo $PATH`). This way you could also add all the other tools to have them at you convenience on your path without having to type the full path to access them.
  ```
  mkdir ~/bin
  ln -s /home/bioinformatics/software/bowtie/bowtie-1.2.1.1/bowtie ~/bin/.
  ln -s /home/bioinformatics/software/bowtie/bowtie-1.2.1.1/bowtie-build ~/bin/.
  ln -s /home/bioinformatics/software/bowtie/bowtie-1.2.1.1/bowtie-inspect ~/bin/.
  ```

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data.
> - Create an new `alignment_tophat` directory
> - Create a new `job.sh` to run tophat on the small fastQ file `SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq`
> - Send job to the cluster
> - Check the output while running and your job using `sacct` on the cluster
> - Update your `README.txt` file with what you've done while job is running
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Viewing alignment data

- Aligned sequence formats: SAM/BAM/CRAM
- SAM format stands for Sequence Alignment Map format and it is the standard format for aligned sequence data which is recognised by the majority of software and browsers

Explore alignment data file format from [SAM/BAM and related specifications](http://samtools.github.io/hts-specs/) and the [SAM Format Specification](http://samtools.github.io/hts-specs/SAMv1.pdf).

To view the header of your aligned reads, you can use `samtools view -H my_file.bam` or to view some aligned reads use `samtools view my_file.bam | tail -10`.

- See [samtools](http://www.htslib.org/) document for more information.
- Our installed version is located here `/home/bioinformatics/software/samtools/samtools-1.6/bin/samtools`.
- To add it onto your path type `ln -s /home/bioinformatics/software/samtools/samtools-1.6/bin/samtools ~/bin/.`

An extract of a BAM read:
```
NS500222:320:HHMJ3BGX3:1:11110:23779:10874	16	8	79810356	60	75M	*	0	0	CAAGGATGCAGCTGTTAGCCCCAATATTTTTTTTATCTTCTCTTGGCACTTTCCTTCACATCTCCTCAGTTGTGC	EEE<EEEEEEAEAEEEEEEEEEAEEEEEEEEEAEEEEAAEAEAEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA	NM:i:0	MD:Z:75	AS:i:75	XS:i:21	RG:Z:1
```
Each aligned reads have:
- column 1: a read name
- column 2: a read alignment flag, go to [Decoding SAM flags](https://broadinstitute.github.io/picard/explain-flags.html)
- column 3: a chromosome to which the read aligns
- column 4: a position in chromosome to which the read aligns
- column 5: a mapping quality
- column 6: alignment information (Cigar string e.g. 100M stands for continuous match of 100 bases)
- column 10: the sequence
- column 11: encoded sequencing quality
- last columns are user defined flags including read group reference defined when running the alignment step

Before being able to visualised your aligned reads, we will need to sort them and create an index, we can use `samtools sort` and `samtools index`. Then download [IGV](http://software.broadinstitute.org/software/igv/) to visualised your BAM files.

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data and into your `alignment` directory.
> - Create a new `job_samtools.sh` to run `samtools sort` and `samtools index` on the small bam file generated previously `SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam`
> - Send job to the cluster
> - Check the output while running and your job using `sacct` on the cluster
> - Update your `README.txt` file with what you've done while job is running
> - Install [IGV](http://software.broadinstitute.org/software/igv/) on you local computer
> - Mount your scratch space onto your local computer and visualise your aligned reads in IGV
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Sort and mark duplicates using Picard tools

- See [Picard tools](http://broadinstitute.github.io/picard/) documentation.
- Our installed version is located here `/home/bioinformatics/software/picard/picard-2.14.0/picard.jar`, you need `java` installed
- Running Picard tools:
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

> :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
>
> - Navigate to your project data and into your `alignment` directory.
> - Create a new `job_picard.sh` to run `SortSam` and `MarkDuplicates` on the small bam file generated previously `SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam`
> - Send job to the cluster
> - Check the output while running and your job using `sacct` on the cluster
> - Update your `README.txt` file with what you've done while job is running
>
> :tada: Congratulations! :thumbsup: You did it! :wink:


## Alignment quality metrics using Picard tools

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

  > :computer: **EXERCISE** Go to your Terminal window, or open a new one and log in onto the cluster head node.
  >
  > - Navigate to your project data and into your `alignment` directory.
  > - Create a new `job_picardmetrics.sh` to run `CollectAlignmentSummaryMetrics` and `CollectInsertSizeMetrics` on the small bam file generated previously `SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam`
  > - Send job to the cluster
  > - Check the output while running and your job using `sacct` on the cluster
  > - Update your `README.txt` file with what you've done while job is running
  >
  > :tada: Congratulations! :thumbsup: You did it! :wink:


## Reference materials

- [CRUK Bioinformatics Autumn School 2017: Functional Genomics](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/)
  - [Short Read Alignment](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/Introduction/SS_DB/Materials/Lectures/Lecture3_ShortRead_Alignment_SS.pdf)
