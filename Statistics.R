require(lattice)
require(ggplot2)

#Read meta data and get grouped samples
Samples <- read.csv("./gse_matrix2.csv", header = TRUE)
Samples <- Samples[,-c(1:40)]
Samples <- Samples[,-c(2,3)]

#C1 = Cluster 1, C2 = Cluster 2
#Each sample's cluster according to Kmeans
Kmeans<-c("C1","C1","C1","C1","C1","C1","C2",
          "C1","C1","C2","C1","C2","C2","C2",
          "C2","C2","C1","C1","C1","C2","C2",
          "C2","C1","C1","C2","C1","C1","C2",
          "C1","C1","C2","C1","C1","C1","C1",
          "C1","C1","C1","C1","C1","C1","C1",
          "C1","C2","C2","C2","C1","C2","C1",
          "C1","C1","C1","C2")

#Chi-square of Samples vs. Kmeans
chisq.test(Samples,Kmeans, correct=FALSE)

#Each sample's cluster according to Hclust
Hclust<-c("C1","C1","C1","C1","C1","C1","C2",
          "C1","C1","C2","C1","C2","C2","C2",
          "C2","C2","C1","C1","C1","C2","C2",
          "C1","C1","C1","C2","C1","C1","C1",
          "C1","C1","C2","C1","C1","C1","C1",
          "C1","C1","C1","C1","C1","C1","C1",
          "C1","C2","C2","C1","C1","C1","C1",
          "C1","C1","C1","C1")

#Chi-square of Samples vs. Hclust
chisq.test(Samples,Hclust, correct=FALSE)

#Each sample's cluster according to PAM
PAM<-c("C1","C1","C1","C1","C1","C1","C2",
       "C1","C1","C2","C1","C2","C2","C2",
       "C2","C2","C1","C1","C1","C2","C2",
       "C2","C1","C1","C2","C1","C1","C1",
       "C1","C1","C2","C1","C1","C1","C1",
       "C1","C1","C1","C1","C1","C1","C1",
       "C1","C2","C2","C2","C1","C2","C1",
       "C1","C1","C1","C2")

#Chi-square of Samples vs. PAM
chisq.test(Samples,PAM, correct=FALSE)

#Chi-square of Kmeans vs. Hclust
chisq.test(Kmeans,Hclust, correct=FALSE)

#Chi-square of Kmeans vs. PAM
chisq.test(Kmeans,PAM, correct=FALSE)

#Chi-square of Hclust vs. PAM
chisq.test(Hclust,PAM, correct=FALSE)

#P values fro chi-squared tests
p_values <-c(0.159,0.04875,0.06262,7.152e-09,3.036e-12,1.54e-09)
#P values adjusted according to Benjamini Hochberg method
p_values_adjusted <- p.adjust(p_values,method="fdr")

Tests <- c("Samples/Kmeans", "Samples/Hclust", "Samples/PAM",
           "Kmeans/Hclust", "Kmeans/PAM", "Hclust/PAM")

#Statistical test results table
stat_test_results <- data.frame(Tests,p_values,p_values_adjusted)

#Enrichment plot
pairs(stat_test_results[2:3], pch = 21)
