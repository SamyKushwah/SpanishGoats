library("readxl")
library("dplyr")
library("ggplot2")
library("matrixStats")

#read in all 4 tables

my_data3 <- read_excel("./FPKM_Data/Data3.xlsx")
#View(my_data3)

my_data4 <- read_excel("./FPKM_Data/Data4.xlsx")
#View(my_data4)

my_data5 <- read_excel("./FPKM_Data/Data5.xlsx")
#View(my_data5)

my_data6 <- read_excel("./FPKM_Data/Data6.xlsx")
#View(my_data6)

#merge B and C tables in their own tables

table_b <- union_all(my_data6,my_data4)
#View(table_c)

table_c <- union_all(my_data5,my_data3)
#View(table_c)

#merge B and C into a single expression matrix
merge_all <- union_all(table_b, table_c)
#View(merge_all)

expression_matrix <- merge_all[!duplicated(merge_all$gene_id), ]
View(expression_matrix)

#replace all n/a values with 0
expression_matrix[is.na(expression_matrix)] = 0
#finds the max values in each row(gene) and adds it to the end of expression_matrix
expression_matrix$row_maximum = rowMaxs(as.matrix(expression_matrix[,c(-1)]))
expression_matrix

#creates a data frame of only the max values
row_range <- expression_matrix[ , ncol(expression_matrix), drop = FALSE]

#plotting the density 
ggplot(row_range, aes(log2(x = row_maximum))) + geom_density(alpha = .2) +ggtitle("Density Plot")
