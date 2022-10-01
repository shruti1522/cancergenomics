func = function(x){
  m =read.csv("C:/Users/grewa/Documents/Python/can_gen1.csv",sep=",",header=T,row.names = 1)
  x <- as.numeric(unlist(m))
  mat <- as.matrix(x)
  for (i in 1:ncol(mat)) {   
    mat[,i] = (mat[,i]/sum(mat[,i]))*1000000
    print(head(mat))
    mat[,i]= log2(mat[,i] +1)
    logfc=log2(mat+1)
    
  }
  return((Heatmap(mat)[1:10]))
}
data = func(x)
data

library(ComplexHeatmap)
library(matrixStats)
data =read.csv("C:/Users/grewa/Documents/Python/can_gen1.csv",sep=",",header=T,row.names = 1)
zscore =function(data){
  x <- as.numeric(data)
  mat <- as.matrix(data)
  for (i in 1:ncol(mat)){
    z_score = (mat- rowMeans(mat))/rowSds(as.matrix(mat))[row(mat)]
    
  }
  z_score[is.na(z_score)]=0
  zcs = as.matrix(z_score)
  return(Heatmap(zcs))


}
zscore(data)



