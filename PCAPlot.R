library("readxl")
library("dplyr")
library(BiocManager)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(tibble)
library(xlsx)

#loading the created expression matrix
expression_matrix <- read.csv("./expressionMatrix.csv", header = TRUE)
View(expression_matrix)

#creating count_data to be used in dds plot object
count_data = expression_matrix[, 2:55]
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
meta_matrix <- meta_matrix %>%
  dplyr::mutate(
    treatment.ch1 = factor(treatment.ch1, levels = c("non transported goats",
                                                         "180 minutes of transportation",
                                                         "30 minutes of transportation"))
  )
levels(meta_matrix$treatment.ch1)
#plot object
vds<-vst(dds, blind=FALSE)

#plotting PCA
plotPCA(vds, intgroup=c("treatment.ch1")) + ggtitle("PCA Plot")

#=======
#=======

#T-SNE Plot
library(M3C)
tsne(count_data, label=as.factor(dds$treatment.ch1))

#=======
#=======

#volcano plot
#doing the differential analysis and storing it
dseq_obj <- DESeq(dds)
dseq_results <- results(dseq_obj)
head(dseq_results)

#make dataset into data frame
deseq_df <- dseq_results %>%
  as.data.frame()
head(deseq_df)

#creating the volcano plot
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = rownames(deseq_df),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
#plotting the volcano plot
volcano_plot

#=======
#=======
library("readxl")
library("dplyr")
library(BiocManager)
library(DESeq2)
library(tidyverse)
library(ggplot2)

#loading the created expression matrix
expression_matrix <- read.csv("./expressionMatrix.csv", header = TRUE)
View(expression_matrix)

#creating count_data to be used in dds plot object
count_data = expression_matrix[, 2:55]
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


#Table of differentially expressed genes 
gene_names <- expression_matrix[,1, drop = FALSE]
diff_exp <- cbind(gene_names, deseq_df)

#find genes that have a p < .05
significant_genes <- subset(diff_exp, pvalue < .05)

#export the differentially expressed genes into their own table
write.xlsx(significant_genes, "significant_genes.xlsx")


#Method 1 gProfiler2
library(gprofiler2)
#enrichement analysis using gprofiler

#get all the gene_ids of the statistically signigicant genes
gene_list_names <- significant_genes$gene_id

#use gostres function to do enrichment analysis 
gostres <- gost(query = gene_list_names, 
                organism = "chircus", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

#use gostplot to see anaylsis
gostplot(gostres, capped = TRUE, interactive = TRUE)









#this method doesn't work
#Method 2: clustProfiler
# http://guangchuangyu.github.io/2015/02/kegg-enrichment-analysis-with-latest-online-data-using-clusterprofiler/
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# https://www.genome.jp/kegg/catalog/org_list.html

#library(clusterProfiler)
#library(DOSE)
#library(MOMA)

#gene_list_names <- significant_genes$gene_id

#significant_genes_entrez <- mapHugo(gene_list_names)
#view(significant_genes_entrez)

#need to change to entrez gene ID 
#https://rdrr.io/bioc/MOMA/man/mapHugo.html
#our gene IDs do not match what they expect

#chxKEGG = enrichKEGG(significant_genes_entrez, organism="chx")

#kk <- enrichKEGG(significant_genes_entrez, organism="chx", pvalueCutoff=0.05, pAdjustMethod="BH", 
#                 qvalueCutoff=0.1)

#head(summary(kk))

#library(enrichplot)

#organism = "human.db"
#BiocManager::install(organism, character.only = TRUE)
#library(organism, character.only = TRUE)

#gse <- gseGO(geneList=gene_list, 
#             ont ="ALL", 
#             keyType = "ENSEMBL", 
#             nPerm = 10000, 
#             minGSSize = 3, 
#             maxGSSize = 800, 
#             pvalueCutoff = 0.05, 
#             verbose = TRUE, 
#             OrgDb = organism, 
#             pAdjustMethod = "none")

#require(DOSE)
#dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)




