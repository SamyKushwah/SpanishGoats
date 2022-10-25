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
library(ggalluvial)

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

k_means <- function(matrix, k_chosen) {
  #Method 1: k means

  km2 <- kmeans(matrix, centers = k_chosen, nstart = 100)
  fviz_nbclust(matrix, kmeans, method ="wss")
  fviz_cluster(kmeans(matrix, centers = k_chosen, nstart = 100), data = matrix)
}

h_clust <- function(matrix) {
  #Method 2: hclust
  
  matrix <- na.omit(matrix)
  matrix <- scale(matrix)
  
  d <- dist(matrix, method = "euclidean")
  
  hc1 <- hclust(d, method = "complete" )
  
  plot(hc1)
  
}

PAM <- function(matrix, k_chosen) {
  #Method 3: PAM Clustering
  
  fviz_nbclust(matrix, pam, method ="silhouette")+theme_minimal()
  pamResult <-pam(matrix, k = k_chosen)
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


#Run with 5000 genes all three clusters methods
#get matrix with only 5000 genes
matrix_5000 <- subset_matrix(5000, diff_exp)

#k_means function, enter matrix and k value
k_means(matrix_5000, 2)

#h_clust function
h_clust(matrix_5000)

#PAM function
PAM(matrix_5000,2)


#Run with 10000 genes all three clusters methods
#get matrix with only 5000 genes
matrix_10000 <- subset_matrix(10000, diff_exp)

#k_means function, enter matrix and k value
k_means(matrix_10000, 2)

#h_clust function
h_clust(matrix_10000)

#PAM function
PAM(matrix_10000,2)


#Run with 1000 genes all three clusters methods
#get matrix with only 5000 genes
matrix_1000 <- subset_matrix(1000, diff_exp)

#k_means function, enter matrix and k value
k_means(matrix_1000, 2)

#h_clust function
h_clust(matrix_1000)

#PAM function
PAM(matrix_1000, 2)


#Run with 100 genes all three clusters methods
#get matrix with only 5000 genes
matrix_100 <- subset_matrix(100, diff_exp)

#k_means function, enter matrix and k value
k_means(matrix_100, 2)

#h_clust function
h_clust(matrix_100)

#PAM function
PAM(matrix_100, 2)


#Run with 10 genes all three clusters methods
#get matrix with only 5000 genes
matrix_10 <- subset_matrix(10, diff_exp)

#k_means function, enter matrix and k value
k_means(matrix_10, 2)

#h_clust function
h_clust(matrix_10)

#PAM function
PAM(matrix_10, 2)




#Alluvial Diagram: 10, 100, 1000, 10000 genes vs clustering methods? 
#amount of clusters? 

cluster_data <- data.frame (k_value = c("2", "2", "2","2", "2", "2","2", "2", "2","2", "2", "2"), 
              clust_method = c("k means","k means","k means","k means","hclust","hclust","hclust",
                               "hclust","PAM Clustering", "PAM Clustering","PAM Clustering", 
                               "PAM Clustering"), num_genes = c("10", "100", "1000", "10000", "10", "100"
                                                                , "1000", "10000",
                                                                "10", "100", "1000", "10000")
)

is_alluvia_form(as.data.frame(cluster_data), axes = 1:3, silent = TRUE)

ggplot(as.data.frame(cluster_data),
       aes(axis1 = num_genes, axis2 = clust_method)) +
  geom_alluvium(aes(fill = k_value), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Number of Genes", "Cluster Method"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Number of Clusters Found Per Cluster Method and Number of Genes")


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

