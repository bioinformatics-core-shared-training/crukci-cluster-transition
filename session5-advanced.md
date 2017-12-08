# Session 5: Advanced analysis and usage of the cluster


## Analysis of RNA-seq data

- Mapping reads on genome using TopHat or STAR
- Transcripts identification and counting using Cufflinks
- Normalisation implemented in edgeR
- Differential expression using DESeq2
- Functional annotation using Blast2GO

### Counting

Once our reads have been aligned against the genome, we need to summarise the information across genes or exons. In the BAM file, there is a chromosomal location for every read that mapped uniquely. We can determine if the region each read is aligned to corresponds to a particular gene or exon and then summarise across the entire BAM file to get total read counts for each gene or exon.

We will use `featureCounts` programme from the [subRead package](http://subread.sourceforge.net/) to do the counting. In addition to the BAM files, we also need to provide `featureCounts` with an annotation file. Usually this will be a GTF/GFF file corresponding to the genome assembly used (a description of the GTF format can be found at [UCSC website](http://genome.ucsc.edu/FAQ/FAQformat.html#format4)). `featureCounts` can also use a simpler annotation format called SAF, this is particularly useful for defining custom/novel features that you wish to count against.

For RNAseq we most commonly wish to count reads aligning to exons, and then to summarise at the gene level. Lets have a quick look at the top of a GTF file so we can see what data it contains and what feature type and attribute type mean:
```
head /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.84.gtf
```

- See `featureCounts` from [subRead package](http://subread.sourceforge.net/) documentation.
- Our installed version is located here `/home/bioinformatics/software/subread/subread-1.5.3/bin/featureCounts`
- Add this tool onto your path `ln -s /home/bioinformatics/software/subread/subread-1.5.3/bin/featureCounts ~/bin/.`
- Running it:
  ```
  featureCounts \
      --primary \
      -C \
      -t exon \
      -g gene_id \
      -a Homo_sapiens.GRCh38.84.gtf \
      -o SLX-14572.i706_i517.featureCounts \
      SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam
  ```
  - `--primary` only count primary alignment
  - `-C` do not count reads where the pairs are mapped to different chromosomes
  - `-t exon` the feature type to count reads against, in this case exons
  - `-g gene_id` the attribute type to summarise counts by, in this case the gene ID

### Filtering out low expressed genes

- See `cpm` function from [edgeR package](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
- To install this package, open RStudio
  ```
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
  ```
- Load it
  ```
  load(limma)
  load(edgeR)
  ```
- Use it:
  ```
  countdata <- read.delim("/Users/pajon01/mnt/scratchb/SLX-14572/featurecounts/SLX-14572.i706_i517.featureCounts", stringsAsFactors = FALSE, comment.char = '#')
  cpmdata <- cpm(countdata)
  head(cpmdata)
  ```

### Differential expression with edgeR

### Adding annotations


## Analysis of ChIP-seq data

### Peak Calling using MACS2

### Peak Annotation

### Differential Binding


## Reference materials

- [CRUK Bioinformatics Autumn School 2017: Functional Genomics](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/)
  - [Introduction to RNA-seq](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/RNASeq/slides/rnaSeq_Sept2017.pdf)
  - [Counting](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/RNASeq/count.nb.html)
  - [Differential Expression Practical](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/RNASeq/slides/LinearModels.pdf)
