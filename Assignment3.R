#Load into R the expression data and matching metadata from GEO that you 
#processed in Assignment 2. 

library("readxl")
library(readr)
library(dbplyr)

expression_matrix <- read_tsv("./GoatDataModified.tsv")
meta_matrix <-read.csv("./gse_matrix2.csv", header = TRUE)
diff_exp <- read_xlsx("./deseq_df.xlsx")

#Get the table of gene names and their expression
gene_names <- expression_matrix[,1, drop = FALSE]
diff_exp <- cbind(gene_names, deseq_df)

#Subset to the 5,000 most variable genes
#Sort by increasing p value
diff_exp <- arrange(diff_exp,pvalue)
#Subset the first 5000 rows
diff_exp <-head(diff_exp, 5000)

