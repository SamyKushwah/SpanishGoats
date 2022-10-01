library("readxl")
library("dplyr")
library(BiocManager)
library(DESeq2)
library(tidyverse)
library(ggplot2)

#loading the created expression matrix
expressionMatrix <- read.csv("./expressionMatrix.csv", header = TRUE)
View(expressionMatrix)

#creating count_data to be used in dds plot object
count_data = expressionMatrix[, 2:55]
#making sure there are no NULL values in the matrix
count_data[is.na(count_data)] = 0
#scaling matrix to make sure there are no 0s by adding 1 to each value
i1 <-which(sapply(count_data, is.numeric))
count_data[i1] <- count_data[i1] +1
View(count_data)

#loading in the metadata 
meta_matrix <-read.csv("./gse_matrix.csv", header = TRUE)
View(meta_matrix)
#removing column for gene_id to make sure col and rows of the 2 matrices match up
col_data = meta_matrix %>% remove_rownames %>% column_to_rownames(var="name")

#making the colnames and rownames of both matrices in the same order
idx = match(rownames(col_data), colnames(count_data))
count_data = count_data[ , idx]

View(count_data)
View(col_data)

#creating an object for the pca plot
dds <- DESeqDataSetFromMatrix(countData = round(count_data),
                              colData = col_data,
                              design = ~ treatment.ch1)

#making sure that the rows and cols are in the same order
all(colnames(count_data) %in% rownames(col_data))

dds

#setting the factor level for comparison from refernce 
dds$treatment.ch1 <- relevel(dds$treatment.ch1, ref = "non transported goats")
dds <- DESeq(dds)

#plot object
vds<-vst(dds, blind=FALSE)

#plotting PCA
plotPCA(vds, intgroup=c("treatment.ch1"))

#T-SNE Plot
library(M3C)
tsne(count_data, label=as.factor(dds$treatment.ch1))

