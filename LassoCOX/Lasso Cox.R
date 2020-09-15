library("glmnet")
library("survival")
rt=read.table("inputfile.txt",header=T,sep="\t",row.names=1)            
rt$futime=rt$futime/365
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit<-glmnet(x, y, family="cox")
plot(fit,label=TRUE)
cvfit=cv.glmnet(x,y,nflod=10, family="cox")
plot(cvfit)
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)