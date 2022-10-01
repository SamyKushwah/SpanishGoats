library("readxl")
library("dplyr")
library("ggplot2")
#read in all 6 tables
#my_data1 <- read_excel("./FPKM_Data/Data1.xlsx")
#View(my_data1)

#my_data2 <- read_excel("./FPKM_Data/Data2.xlsx")
#View(my_data2)

my_data3 <- read_excel("./FPKM_Data/Data3.xlsx")
View(my_data3)

my_data4 <- read_excel("./FPKM_Data/Data4.xlsx")
View(my_data4)

my_data5 <- read_excel("./FPKM_Data/Data5.xlsx")
View(my_data5)

my_data6 <- read_excel("./FPKM_Data/Data6.xlsx")
View(my_data6)

#merge all of the A, B and C tables in their own 3 tables
#tableA <- union_all(my_data1,my_data2)
#View(tableA)

table_b <- union_all(my_data6,my_data4)
View(table_c)

table_c <- union_all(my_data5,my_data3)
View(table_c)

#issues: the three tables don't seem to have any genes in common 
#not sure what to do next
#mergeAB <- union_all(tableA, tableB)
#View(mergeAB)

merge_all <- union_all(table_b, table_c)
View(merge_all)

expression_matrix <- merge_all[!duplicated(merge_all$gene_id), ]
View(expression_matrix)

density_plot_data <- read_excel("./DensityPlotData.xlsx")
View(density_plot_data)
ggplot(density_plot_data, aes(log2(x = x), fill = Sample)) + geom_density(alpha = .2)
