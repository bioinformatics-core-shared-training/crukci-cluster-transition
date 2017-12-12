library(edgeR)
# global setting for data frame to not auto-convert to factors automatically
options(stringsAsFactors = FALSE)

# read sample information
samplesheet <- read.csv("data/samplesheet_RNAseq.csv")
View(samplesheet)
# read count data (output of featureCounts tool)
countdata <- read.csv("data/SLX-12763.AllSamples.featureCounts", comment.char = "#")
View(countdata)

# data manipulations
rownames(countdata) <- countdata$Geneid
countdata <- countdata[,-(1:6)]
samplesheet$CountTableNames <- gsub("[/-]", ".", samplesheet$BamFile)
colnames(countdata) <- samplesheet$SampleName[match(colnames(countdata), samplesheet$CountTableNames)]
View(countdata)

# Filtering out low expressed genes
cpmdata <- cpm(countdata)
head(cpmdata)

# Identify genes with at least 0.5 cpm in at least 2 samples
thresh <- cpmdata > 0.5
keep <- rowSums(thresh) >= 2

# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]

# Convert to an edgeR object
dgeObj <- DGEList(counts.keep)

# Perform TMM normalisation
dgeObj <- calcNormFactors(dgeObj)

# Create the design matrix
groups <- samplesheet$Source
design <- model.matrix(~groups)
View(design)

# Estimating the overall dispersion
dgeObj <- estimateCommonDisp(dgeObj)
dgeObj <- estimateGLMTrendedDisp(dgeObj)

# Estimating gene-wise dispersion estimates
dgeObj <- estimateTagwiseDisp(dgeObj)

# Plot the estimated dispersions
plotBCV(dgeObj)

# Fit the linear model
fit <- glmQLFit(dgeObj, design)
names(fit)
head(coef(fit))

# Conduct likelihood ratio test
qlf <- glmQLFTest(fit, coef=2)
de <- decideTestsDGE(qlf)
summary(de)

# Visualise the results using plotSmear
detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(qlf, de.tags=detags)

# Results table
results <- qlf$table
results$ENSEMBL <- rownames(results)
View(results)

# Adjust p-value
?p.adjust
results$Padj <- p.adjust(results$PValue, method="BH")
results <- results[order(results$Padj),]
View(results)

# Load annotation library
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

# Add gene description to the results table
ann <- select(org.Hs.eg.db,keytype="ENSEMBL", keys=rownames(results), columns=c("ENSEMBL","SYMBOL","GENENAME"))
results <- merge(x=results, y=ann, by= "ENSEMBL", all.x=TRUE)
View(results)
