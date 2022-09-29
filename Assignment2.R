library("readxl")
my_data1 <- read_excel("desktop/FPKM Reads/Data1.xlsx")
View(my_data1)

my_data2 <- read_excel("desktop/FPKM Reads/Data2.xlsx")
View(my_data2)

total <- merge(my_data1,my_data2)
View(total)