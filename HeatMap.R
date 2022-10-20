library(ComplexHeatmap)
library(circlize)
library(readxl)


cancer <-read.csv("E:/MSC/msc 3rd sem 2022/Cancer genomics/can_gen1.csv", row.names = 1)
cancer <- matrix(rnorm(100, 0, 5), nrow = 10, ncol = 10)

colnames(cancer) <- paste0("gene", 1:10)
rownames(cancer) <- paste0("species", 1:10)

Heatmap(cancer)
Heatmap(cancer , name = "mat " , column_title ="Cancer")
col_fun  = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))

Heatmap(cancer, col= col_fun, name = "range", column_title ="Cancer" )




