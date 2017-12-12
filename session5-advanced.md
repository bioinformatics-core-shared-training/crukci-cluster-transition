# Session 5: Advanced analysis and usage of the cluster

## Learning Objectives

- Analysis of RNA-seq data including
  - Run Feature counts
  - Read counts data in R
  - Run differential expression analysis
  - Add annotations
- Overview of tools for analysis of ChIP-seq data


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
      -o SLX-14572.FourSamples.featureCounts \
      SLX-14572.i706_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam SLX-14572.i706_i502.HHMJ3BGX3.s_1.r_1.small.fq.bam SLX-14572.i706_i503.HHMJ3BGX3.s_1.r_1.small.fq.bam SLX-14572.i703_i517.HHMJ3BGX3.s_1.r_1.small.fq.bam
  ```
  - `--primary` only count primary alignment
  - `-C` do not count reads where the pairs are mapped to different chromosomes
  - `-t exon` the feature type to count reads against, in this case exons
  - `-g gene_id` the attribute type to summarise counts by, in this case the gene ID

Running featureCounts generates two output file `SLX-14572.FourSamples.featureCounts` `SLX-14572.FourSamples.featureCounts.summary`.
- The summary table reports the numbers of unassigned reads and the reasons why they are not assigned (eg. ambiguity, multi-mapping, secondary alignment, mapping quality, fragment length, chimera, read duplicate, non-junction and so on), in addition to the number of successfully assigned reads for each library.
- The full results table has multiple columns:
  - column 1: the gene identifier
  - columns 2-5: the genes location
  - column 6: the length of the gene
  - column 7-10: the number of reads assigned to the gene

### Reading the count data in R

- Download and install [RStudio](https://www.rstudio.com/)
- [Get cluster file systems mounted on your macOS](session3-cluster-usage.md#get-cluster-file-systems-mounted-on-your-macos)
- Download [GitHub crukci-cluster-transition data](https://github.com/bioinformatics-core-shared-training/crukci-cluster-transition/archive/master.zip)
- Unzip `crukci-cluster-transition-master.zip`

- Set up an RStudio project:
  - Menu *File* > *New Project...*
  - *New Directory* > *Empty Project*
  - **Directory name:** type **.**
  - **Create project as subdirectory of:** click *Browse...*
  - and select the directory `crukci-cluster-transition-master`, click *Open*
  - click *Create Project*

- Learn R with [Half-day introduction to the R language crash course](https://bioinformatics-core-shared-training.github.io/r-crash-course/)

- Install package dependencies for RNAseq analysis in R
  ```R
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
  biocLite("org.Mm.eg.db")
  ```

- Read sample information in R
  ```R
  samplesheet <- read.csv("data/samplesheet.SLX-12345.csv")
  View(samplesheet)
  ````

- Read count data in R
  ```R
  countdata <- read.csv("data/SLX-12345.AllSamples.featureCounts", comment.char = "#")
  View(countdata)
  ```

### Data manipulations

Convert Gene names into row names, remove columns, and rename columns
```R
rownames(countdata) <- countdata$Geneid
countdata <- countdata[,-(1:6)]
samplesheet$CountTableNames <- gsub("-", ".", samplesheet$BamFile)
colnames(countdata) <- samplesheet$SampleName[match(colnames(countdata), samplesheet$CountTableNames)]
```

### Filtering out low expressed genes

Genes with very low counts across all libraries provide little evidence for differential expression and they interfere with some of the statistical approximations that are used later in the pipeline.

- See `cpm` function from [edgeR package](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
- To install this package in RStudio
  ```R
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
  ```
- Load it
  ```R
  library(edgeR)
  ```

- Calculate the Counts Per Million measure
  ```R
  cpmdata <- cpm(countdata)
  head(cpmdata)
  ```
- Identify genes with at least 0.5 cpm in at least 2 samples
  ```R
  thresh <- cpmdata > 0.5
  keep <- rowSums(thresh) >= 2
  ```
- Subset the rows of countdata to keep the more highly expressed genes
  ```R
  counts.keep <- countdata[keep,]
  ```

### Converting counts to DGEList object

- Convert to an edgeR object
  ```R
  dgeObj <- DGEList(counts.keep)
  ```

- Perform TMM normalisation
  ```R
  dgeObj <- calcNormFactors(dgeObj)
  ```

### Create the design matrix

```R
# Create the groups
groups <- samplesheet$Group
# Specify a design matrix
design <- model.matrix(~groups)
```

### Differential expression with edgeR

- Estimating the overall dispersion
  ```R
  dgeObj <- estimateCommonDisp(dgeObj)
  ```
- Estimating gene-wise dispersion estimates
  ```R
  dgeObj <- estimateGLMTrendedDisp(dgeObj)
  dgeObj <- estimateTagwiseDisp(dgeObj)
  ```
- Plot the estimated dispersions
  ```R
  plotBCV(dgeObj)
  ```
- Fit the linear model
  ```R
  fit <- glmFit(dgeObj, design)
  names(fit)
  head(coef(fit))
  ```
- Conduct likelihood ratio test
  ```R
  qlf <- glmQLFTest(fit, coef=2)
  de <- decideTestsDGE(qlf)
  summary(de)
  ```
- Visualise the results using plotSmear
  ```R
  detags <- rownames(dgeObj)[as.logical(de)]
  plotSmear(qlf, de.tags=detags)
  ```

### Output table

- Extract table from `qlf` object
```R
results <- qlf$table
results$ENSEMBL <- rownames(results)
View(results)
```
- Adjust p-value
```R
results$Padj <- p.adjust(results$PValue, method="BH")
results <- results[order(results$Padj),]
View(results)
```

### Volcano plot

```R
deg <- which(results$Padj<=0.05)
plot(results$logFC, -log10(results$Padj), pch=21, col="black", bg="black", xlab="log2(FoldChange)", ylab="-log10(adjusted p-value)")
points(results$logFC[deg], -log10(results$Padj[deg]), pch=21, col="black", bg="red")
```

### Add annotations

There are a number of ways to add annotation, but we will demonstrate how to do this using the `org.Mm.eg.db` package. This package is one of several organism-level packages which are re-built every 6 months. These packages are listed on the [annotation section](http://bioconductor.org/packages/release/BiocViews.html#___AnnotationData) of the Bioconductor, and are installed in the same way as regular Bioconductor packages.

- To install this package in RStudio
  ```R
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Mm.eg.db")
  ```
- Load it
  ```R
  library(org.Mm.eg.db)
  ```
- View columns
  ```R
  columns(org.Mm.eg.db)
  ```
- Add gene description to the results table
  ```R
  ann <- select(org.Hs.eg.db,keytype="ENSEMBL", keys=rownames(results), columns=c("ENSEMBL","SYMBOL","GENENAME"))
  results <- merge(x=results, y=ann, by= "ENSEMBL", all.x=TRUE)
  View(results)
  ```

## Analysis of ChIP-seq data

- Mapping reads on genome using bwa
- Quality Control using ChIPQC
- Peak Calling with MACS2
- Peak Annotation using ChIPseeker
- Motif Analysis using Meme Suite
- Normalisation implemented in DiffBind
- Differential Binding Analysis implemented in DiffBind using DeSeq2

### Peak Calling using MACS2

Once our reads have been aligned against the genome, we need to identify regions of enrichment (peaks). There are a variety of tools
available for calling peaks: SICER, MACS2, EPIC, Enriched Domain Detector (EDD), BayesPeak etc. Here we will use MACS2.

Different  types of ChIP data have differently shaped  peaks. Generally, TF peaks are narrow, whilst epignomic data, such as histone marks, can be narrow, broad, or a mixture of both. It is important to use a peak caller that is appropriate for the peak type being sought. MACS2 has both narrow and broad modes and so is widely applicable.

When calling peaks for a sample it is also necessary to provide an appropriate input (control) sample. This is a negative control that will allow the peak caller to estimate the background signal. Some peak callers will work without an input sample but this is **not** recommended.

- MACS2 can be downloaded [here](https://github.com/taoliu/MACS)
- Our installed version is located here **...** (not yet installed)
- Add this tool onto your path `ln -s /home/bioinformatics/software/... ~/bin/.`
- Running it:
	- Make a directory for the output files `mkdir macs`
	- Run the tool:
```
macs2 callpeak \
    --treatment **ChIPsample.bam** \
    --control **InputSample.bam** \
    --gsize "2671858539" \
    --outdir "macs" \
    --name "ChIPSampleName" `
```
- MACS2 has a large number of options and arguments, the above is a default narrow peak analysis.
- Note the `--gsize` argument - this is the "effective genome size" - this parameter is dependent on the read length and the actual genome size,
please see the MACS documentation for further details.

### Peak Annotation

### Differential Binding


## Reference materials

- [CRUK Bioinformatics Autumn School 2017: Functional Genomics](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/)
  - [Introduction to RNA-seq](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/RNASeq/slides/rnaSeq_Sept2017.pdf)
  - [Counting](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/RNASeq/count.nb.html)
  - [Differential Expression Practical](https://bioinformatics-core-shared-training.github.io/cruk-autumn-school-2017/RNASeq/slides/LinearModels.pdf)
