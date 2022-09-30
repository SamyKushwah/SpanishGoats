library("readxl")
library("dplyr")
#read in all 6 tables
my_data1 <- read_excel("desktop/FPKM Reads/Data1.xlsx")
View(my_data1)

my_data2 <- read_excel("desktop/FPKM Reads/Data2.xlsx")
View(my_data2)

my_data3 <- read_excel("desktop/FPKM Reads/Data3.xlsx")
View(my_data3)

my_data4 <- read_excel("desktop/FPKM Reads/Data4.xlsx")
View(my_data4)

my_data5 <- read_excel("desktop/FPKM Reads/Data5.xlsx")
View(my_data5)

my_data6 <- read_excel("desktop/FPKM Reads/Data6.xlsx")
View(my_data6)

#merge all of the A, B and C tables in their own 3 tables
tableA <- union_all(my_data1,my_data2)
View(tableA)

tableB <- union_all(my_data4,my_data6)
View(tableB)

tableC <- union_all(my_data3,my_data5)
View(tableC)

#issues: the three tables don't seem to have any genes in common 
#not sure what to do next
mergeAB <- union_all(tableA, tableB)
View(mergeAB)

mergeAll <- union_all(mergeAB, tableC)
View(mergeAll)

