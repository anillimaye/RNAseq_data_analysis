# Loading the required libraries

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
# Load the countdata
```
Data <- read.table(file.choose(), header=TRUE, sep="\t")
head(Data)
View(Data)
```
# Removing extra Geneid columns
```
Data <- Data[,-c(3,5,7,9,11,13)]
```
#renaming the column names that are more friendly
names(Data) <- c("ensgene","C1","C2","C3","EL1","EL2","EL3","EL4")

#now we will add rownames to the data, which will actually be the ensgene column
rownames(Data) <- Data$ensgene

#Now we dont want the ensgene column
Data <- Data[, -1]

#We will see summary of the data which shows that large number of entries are zeros
#in each column
summary(Data)

#One can see how many total reads are there in each sample
#and visualize them as bar plots
colSums(Data)
barplot(colSums(Data), las=3)

#We look at the distribution of counts in each column(sample)
hist(Data$C1, breaks=100)

#log2 transformation of the data
logData <- log2(1+Data)

#Now you can see the change in the histogram
hist(logData$C1, br=100)

#at this point one can look at similarity in the logcountdata in any 
#pair of samples; for example

plot(logData$C1, logData$C2)

#We will construct the metadata of the experiment
Line <- c("MCF7","MCF7","MCF7","MCF7","MCF7","MCF7","MCF7")
Treatment <- c("Control", "Control", "Control", "Enternolactone", "Enternolactone", "Enternolactone", "Enternolactone")
Time <- c("24", "24", "24", "24", "24", "24", "24")
Sample <- names(Data)

colData <- as.data.frame(cbind(Sample, Time, Line, Treatment))

#Now we make the DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=Data, colData=colData, design= ~Treatment)
dds <- DESeq(dds)
dim(dds)

#we remove those ensgenes which have very low counts in all samples
dds <- dds[rowSums(counts(dds))>7, ]
dim(dds)

#now we look at sizeFactors
sizeFactors(dds)

#regularized log transformation of the data
rld <- rlog(dds)

#PCAplot
plotPCA(rld, intgroup="Treatment")
plotPCA(rld, intgroup=c("Sample", "Treatment"))

#generating a heatmap to find out how close are the samples to each other.
pheatmap(cor(assay(rld)))

#this part will determine the differentially expressed genes

#to see how many different kinds of comparisons are possible
resultsNames(dds)
res <- results(dds)
res <- results(dds, alpha=0.05)
head(res)
summary(res)
plotMA(res, ylim=c(-8,8))

#apply lfcshrink
res_lfcshrink=lfcShrink(dds,contrast = c("Treatment","Enternolactone","Control"),type = 'normal',res = res)
plotMA(res_lfcshrink, ylim=c(-8,8))

#subsetting the lfcshrinked results to keep only those genes where padj < 0.05
res_lfcshrink_sig=subset(res_lfcshrink,padj<0.05)

#saving the lfcshrinked results into a new object name
sig_res <- res_lfcshrink_sig

#getting gene symbols against the ensembl ids#first we change the rownams to column for the object sig_res
sig_res <- data.frame(sig_res) %>% rownames_to_column(var="ensgene")
#then we make a list of ensembl ids
ensembl_ids <- sig_res$ensgene

listEnsembl()
ensembl <- useEnsembl(biomart="genes")
datasets <- listDatasets(ensembl)
#connect ensembl
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

#connect ensembl
ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#get symbol against each ensembl id

pair <- getBM(attributes=c("ensembl_gene_id","external_gene_name"), filters="ensembl_gene_id", values=ensembl_ids, mart=ensembl.con)

# changing the column headers of pair to "ensgene" and "symbol"
names(pair) <- c("ensgene","symbol")

## Merging the two objects sig_res and pair using the merge command
```
final <- merge(sig_res, pair, by="ensgene")
```

