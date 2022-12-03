### Loading the required libraries

First load all the libraries that will be needed for the processing and analysis of data.
```
library(DESeq2)
library(dplyr)
library(ggplot2)
library(gplots)
library(annotables)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)
library(data.table)
```
### Load the countdata

Upon executing the first line of the code block below, you will select the file from your system that contains the raw count data. 
```
Data <- read.table(file.choose(), header=TRUE, sep="\t")
head(Data)
View(Data)
```
### Removing extra Geneid columns
```
Data <- Data[,-c(3,5,7,9,11,13)]
```
### Renaming the column names that are more friendly
```
names(Data) <- c("ensgene","C1","C2","C3","EL1","EL2","EL3","EL4")
```

### Add rownames to the data, which will actually be the ensgene column
```
rownames(Data) <- Data$ensgene
```
### Removal of ensgene columns since rownames have been included in the previous step
```
Data <- Data[, -1]
```
### Display summary of the data
```
summary(Data)
```

### Determining total reads in each of the samples and visualizing them as barplots
```
colSums(Data)
barplot(colSums(Data), las=3)
```
### To check the distribution of counts for each sample
hist(Data$C1, breaks=100)

### log2 transformation of the data
```
logData <- log2(1+Data)
```
### Distribution of log transformed counts in each sample
```
hist(logData$C1, br=100)
```
### Examining a pair of samples for similarity in the logcountdata in any 
```
plot(logData$C1, logData$C2)
```

### Constructing the metadata
```
Line <- c("MCF7","MCF7","MCF7","MCF7","MCF7","MCF7","MCF7")
Treatment <- c("Control", "Control", "Control", "Enternolactone", "Enternolactone", "Enternolactone", "Enternolactone")
Time <- c("24", "24", "24", "24", "24", "24", "24")
Sample <- names(Data)

colData <- as.data.frame(cbind(Sample, Time, Line, Treatment))
```
### Now we make the DESeq2 object
```
dds <- DESeqDataSetFromMatrix(countData=Data, colData=colData, design= ~Treatment)
dds <- DESeq(dds)
dim(dds)
```
### We remove those ensgenes which have very low counts in all samples
```
dds <- dds[rowSums(counts(dds))>7, ]
dim(dds)
```
### Looking at sizeFactors
```
sizeFactors(dds)
```
### Regularized log transformation of the data
```
rld <- rlog(dds)
```
### PCAplot
```
plotPCA(rld, intgroup="Treatment")
plotPCA(rld, intgroup=c("Sample", "Treatment"))
```
### Generating a heatmap to find out how close are the samples to each other.
```
pheatmap(cor(assay(rld)))
```
## This part will determine the differentially expressed genes

### To see how many different kinds of comparisons are possible
```
resultsNames(dds)
res <- results(dds)
res <- results(dds, alpha=0.05)
head(res)
summary(res)
plotMA(res, ylim=c(-8,8))
```
### Apply lfcshrink
```
res_lfcshrink=lfcShrink(dds,contrast = c("Treatment","Enternolactone","Control"),type = 'normal',res = res)
plotMA(res_lfcshrink, ylim=c(-8,8))
```
### Subsetting the lfcshrinked results to keep only those genes where padj < 0.05
```
res_lfcshrink_sig=subset(res_lfcshrink,padj<0.05)
sig_res <- res_lfcshrink_sig
```
### Getting gene symbols against the ensembl ids#first we change the rownams to column for the object sig_res
```
sig_res <- data.frame(sig_res) %>% rownames_to_column(var="ensgene")
# then we make a list of ensembl ids
ensembl_ids <- sig_res$ensgene

listEnsembl()
ensembl <- useEnsembl(biomart="genes")
datasets <- listDatasets(ensembl)
# connect ensembl
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

#connect ensembl
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#get symbol against each ensembl id
pair <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), filters="ensembl_gene_id", values=ensembl_ids, mart=ensembl.con)

#changing the column headers of pair to "ensgene" and "symbol"
names(pair) <- c("ensgene","symbol")
```
## Merging the two objects sig_res and pair using the merge command
```
final <- merge(sig_res, pair, by="ensgene")
```

