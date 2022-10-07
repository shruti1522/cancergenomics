p=read.csv("E:/MSC/msc 3rd sem 2022/Cancer genomics/can_gen.xlsx",sep=",",header=T,row.names = 1)
p



function1=function(p){
  matrixc=p
  m <- as.numeric(unlist(matrixc))
  for (i in 1:ncol(x)) {
    matrixc[,i] = (x[,i]/sum(x[,i]))*1000000
    print(head(matrixc))
    matrixc[,i]= log2(matrixc[,i] +1)
    logfc=log2(matrixc+1)
  }
  mat=matrix(NA,ncol=4,nrow = nrow(logfc))
  rownames(mat)=rownames(logfc)
  colnames(mat)=c('meanTumor','meanControl','pvalue','log2FC')
}
}
for(i in 1:nrow(mat)){
  vector1 = as.numeric(mat[i, 1:2])
  
  
  vector2 = as.numeric(mat[i, 2:4])
  
  res=t.test(vector1, vector2, paired = F, alternative = "two.sided")
  mat[i,1]=res$estimate[[1]]
  mat[i,2]=res$estimate[[2]]
  mat[i,3]=res$p.value
  mat[i,4]=mat[i,1]-mat[i,2]
  
}
}
mat=as.data.frame(mat)
num=which(is.nan(mat$pvalue))
mat[num,'pvalue']=1

library(EnhancedVolcano)
EnhancedVolcano(mat,lab = rownames(mat),x = 'log2FC' ,y ='pvalue')