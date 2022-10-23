#Load into R the expression data and matching metadata from GEO that you 
#processed in Assignment 2. 

library("readxl")
library(readr)
library(dplyr)
library(gtools)

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
#Join the genes from diff_exp to a gene matrix
significant_genes <- data.frame(diff_exp$Gene)
colnames(significant_genes)[1] = "Gene"

expression_matrix_5000 <- inner_join(significant_genes,expression_matrix, by = "Gene")





#Method 1:ConsensusClusterPlus
library("ConsensusClusterPlus")

matrix_5000 <- expression_matrix_5000
matrix_5000<-matrix_5000[,-1] # delete column of genes
colnames(matrix_5000)<-NULL 
matrix_5000 <- as.matrix(matrix_5000)

rcc = ConsensusClusterPlus(matrix_5000,maxK=7,reps=1000,pItem=0.8,pFeature=1,
                           title="example",distance="pearson",clusterAlg="hc")
#resICL = calcICL(rcc,title="example")
calcICL(rcc,title="untitled_consensus_cluster")

#Method 2: Gaussian Mixture Models
library(ClusterR)

opt_gmm = Optimal_Clusters_GMM(matrix_5000, max_clusters = 7, criterion = "BIC", 
                               
                               dist_mode = "maha_dist", seed_mode = "random_subset",
                               
                               km_iter = 10, em_iter = 10, var_floor = 1e-10, 
                               
                               plot_data = T)

#Method 3: k means
library(ggplot2)
library(stats)
library(factoextra)

km2 <- kmeans(matrix_5000, centers = 5, nstart = 100)
fviz_nbclust(matrix_5000, kmeans, method ="wss")
fviz_cluster(kmeans(matrix_5000, centers = 5, nstart = 100), data = matrix_5000)

#Method 4: hclust
library(tidyverse)  
library(cluster)    
library(factoextra) 
library(dendextend) 

matrix_5000 <- na.omit(matrix_5000)
matrix_5000 <- scale(matrix_5000)

d <- dist(matrix_5000, method = "euclidean")

hc1 <- hclust(d, method = "single" )

plot(hc1)


#Method 5: PAM
library(cluster)
pamx <- pam(matrix_5000, 7)

