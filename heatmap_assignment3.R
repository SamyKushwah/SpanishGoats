library(ComplexHeatmap)
library(DESeq2)
library(ggplot2)

#arranges matrix_100 in alphabetical/group order
matrix_100 <- matrix_100[c(50, 49, 46, 28, 27, 47, 26, 
                           25, 24, 43, 42, 52, 40, 39, 
                           38, 41, 37, 36, 14, 13, 12, 
                           11, 10, 48, 53, 45, 44, 8, 
                           6, 5, 4, 7, 3, 2, 1, 
                           23, 22, 21, 20, 19, 18, 17,
                           16, 15, 35, 34, 32, 31, 51,
                           33, 30, 29, 9),]

mat <- as.data.frame(matrix_100)

mat = t(mat)

#create vector for sample grouping
li = c("Control")
for (x in 1:17)
  li <- append(li, "Control")
for (x in 1:17)
  li <-append(li, "30 min")
for (x in 1:18)
  li <-append(li, "180 min")

#C1 = Cluster 1, C2 = Cluster 2
#Each sample's cluster according to Kmeans
Kmeans<-c("C1","C1","C1","C1","C1","C1","C2",
          "C1","C1","C2","C1","C1","C2","C1",
          "C1","C1","C1","C1","C1","C2","C2",
          "C2","C1","C1","C2","C1","C1","C2",
          "C1","C2","C1","C1","C1","C1","C2",
          "C2","C1","C2","C2","C2","C2","C2",
          "C1","C1","C1","C1","C1","C1","C1",
          "C1","C1","C2","C2")


#Each sample's cluster according to Hclust
Hclust<-c("C1","C1","C1","C1","C1","C1","C2",
          "C1","C1","C1","C1","C1","C2","C1",
          "C1","C1","C1","C1","C1","C2","C2",
          "C1","C1","C1","C2","C1","C1","C1",
          "C1","C1","C1","C1","C1","C1","C1",
          "C2","C1","C2","C2","C2","C2","C2",
          "C1","C1","C1","C1","C1","C1","C1",
          "C1","C1","C2","C2")


#Each sample's cluster according to PAM
PAM<-c("C1","C1","C1","C1","C1","C1","C2",
       "C1","C1","C1","C1","C1","C2","C1",
       "C1","C1","C1","C1","C1","C2","C2",
       "C2","C1","C1","C2","C1","C1","C2",
       "C1","C2","C1","C1","C1","C1","C2",
       "C2","C1","C2","C2","C2","C2","C2",
       "C1","C1","C1","C1","C1","C1","C1",
       "C1","C1","C2","C2")

#Top annotation
ha = HeatmapAnnotation(Group = li, border = T)

#Bottom annotation
hb = HeatmapAnnotation(Kmeans_Clusters = Kmeans,
                       Hclust_Clusters = Hclust,
                       PAM_Clusters = PAM,
                       border = T)

#Creates Heatmap
Heatmap(mat, 
        cluster_rows = T, cluster_columns = T, 
        column_labels = colnames(mat), 
        row_names_gp = grid::gpar(fontsize = 6),
        name = "Counts", 
        top_annotation = ha, bottom_annotation = hb)
