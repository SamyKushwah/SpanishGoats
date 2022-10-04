library(devtools)
library("readxl")
library("dplyr")
library(BiocManager)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(tibble)
#BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
#getting the names of significant genes
test <- significant_genes[order(significant_genes$pvalue),]
test2 <- data.frame(test$gene_id)
test3 <- head(test2, 70)

#omitting duplicates
test3 <- test3[!duplicated(test3$test.gene_id), ]
test3 <- test3[-c(6)]

test3 <- as.data.frame(test3)


#omitting duplicates from expression matrix
duplicate_text <- expression_matrix[!duplicated(expression_matrix$gene_id), ]
duplicate_text <- duplicate_text[-c(3), ]

#inner join between expression matrix and significant genes
result <- merge (x=test3, 
                 y=duplicate_text, 
                 by.x = "test3",
                 by.y = "gene_id",
                 x.all=FALSE,
                 y.all=FALSE)

colnames(result)[which(names(result) == "test3")] <- "gene_id"

expression_matrix_gene <- result %>% remove_rownames %>% column_to_rownames(var="gene_id")

count_data <- expression_matrix_gene

count_data[is.na(count_data)] = 0

count_data[ , order(names(count_data))]
# list for annotation
x <- c('control','control','control','180 min','control','control','control','control','control','180 min','180 min','180 min','180 min','180 min','180 min','control','180 min','180 min','180 min','180 min','control','180 min','180 min','180 min','control','180 min','180 min','180 min','180 min','control','control','control','control','control','control','control','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min','30 min')
mat = as.matrix(count_data)
ha = HeatmapAnnotation(Group = x, border = T)

#creating heatmap
Heatmap(mat, 
        name = "fpkm reads", 
        column_title = "Heatmap of significant genes determined by p-value",
        top_annotation = ha)










