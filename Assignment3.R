#Load into R the expression data and matching metadata from GEO that you 
#processed in Assignment 2. 

library("readxl")
library(readr)
library(dplyr)
library(gtools)
library(ggplot2)
library(stats)
library(factoextra)
library(tidyverse)  
library(cluster)    
library(factoextra) 
library(dendextend) 

subset_matrix <- function(num, diff_exp) {
  #Subset the first # rows
  diff_exp <-head(diff_exp, num)
  #Join the genes from diff_exp to a gene matrix
  significant_genes <- data.frame(diff_exp$Gene)
  colnames(significant_genes)[1] = "Gene"
  
  expression_matrix <- inner_join(significant_genes,expression_matrix, by = "Gene")
  matrix <- expression_matrix
  matrix<-matrix[,-1] # delete column of genes
  colnames(matrix)<-NULL 
  matrix <- as.matrix(matrix)
  matrix <- t(matrix)
  return(matrix)
}

k_means <- function(matrix) {
  #Method 1: k means
  
  km2 <- kmeans(matrix_5000, centers = 2, nstart = 100)
  fviz_nbclust(matrix_5000, kmeans, method ="wss")
  fviz_cluster(kmeans(matrix_5000, centers = 2, nstart = 100), data = matrix_5000)
}

h_clust <- function(matrix) {
  #Method 2: hclust
  
  matrix_5000 <- na.omit(matrix_5000)
  matrix_5000 <- scale(matrix_5000)
  
  d <- dist(matrix_5000, method = "euclidean")
  
  hc1 <- hclust(d, method = "complete" )
  
  plot(hc1)
  
}

PAM <- function(matrix) {
  #Method 3: PAM Clustering
  
  fviz_nbclust(matrix_5000, pam, method ="silhouette")+theme_minimal()
  pamResult <-pam(matrix_5000, k = 2)
  fviz_cluster(pamResult, 
               palette =c("#007892","#D9455F"),
               ellipse.type ="euclid",
               repel =TRUE,
               ggtheme =theme_minimal())
}

expression_matrix <- read_tsv("./GoatDataModified.tsv")
meta_matrix <-read.csv("./gse_matrix2.csv", header = TRUE)
diff_exp <- read_xlsx("./deseq_df.xlsx")

#Get the table of gene names and their expression
gene_names <- expression_matrix[,1, drop = FALSE]
diff_exp <- cbind(gene_names, deseq_df)

#Sort by increasing p value
diff_exp <- arrange(diff_exp, pvalue)

#get matrix with only 5000 genes
matrix_5000 <- subset_matrix(5000, diff_exp)

#k_means function
k_means(matrix_5000)

#h_clust function
h_clust(matrix_5000)

#PAM function
PAM(matrix_5000)


#Not using these
#Method 4:ConsensusClusterPlus
#library("ConsensusClusterPlus")

#rcc = ConsensusClusterPlus(matrix_5000,maxK=1,reps=50,pItem=0.8,pFeature=1,
#                           title="example",distance="pearson",clusterAlg="hc")

#resICL = calcICL(rcc,title="example")
#calcICL(rcc,title="untitled_consensus_cluster")

#Method 5: Gaussian Mixture Models
#library(ClusterR)

#opt_gmm = Optimal_Clusters_GMM(matrix_5000, max_clusters = 9, criterion = "BIC", 
#                               dist_mode = "maha_dist", seed_mode = "random_subset",
#                               km_iter = 10, em_iter = 10, var_floor = 1e-10, 
#                               plot_data = T)

