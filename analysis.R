library(ComplexHeatmap)
library(matrixStats)
data_file =read.csv("E:/MSC/msc 3rd sem 2022/Cancer genomics/can_gen1.csv",sep=",",header=T,row.names = 1)
data_file

fun =function(data_file){
  matrixc=data_file
  for (i in 1:ncol(data_file)) {
    cpm[,i] = (data_file[,i]/sum(data_file[,i]))*1000000
    print(head(cpm))
    cpm[,i]= log2(cpm[,i] +1)
    logfc=log2(cpm+1)
  }
  d = cpm
  for (i in 1:ncol(data)){
    z_score = (data - rowMeans(d))/rowSds(as.matrix(d))[row(d)]
    
  }
  zscore[is.na(zscore)]=0
  z = as.matrix(zscore)
  
   Heatmap(z)
  return(Heatmap(z[1:10],))
  
}

fun(data_file)