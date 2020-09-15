install.packages('e1071')
library("limma")    
rt=read.table("inputfile.txt",sep="\t",header=T,check.names=F)     
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
data=data[rowMeans(data)>0,]
v <-voom(data, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)   
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)