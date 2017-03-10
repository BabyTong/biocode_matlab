integrate<-function()
{
library("samr")
exp<-read.csv("input/exp.csv",header=FALSE)
exp<-as.matrix(exp)
new_id<-read.csv("input/geneid.csv",header=FALSE)
dimnames(exp)[[1]]<-new_id[[1]]
label<-read.csv("input/label.csv",header=FALSE)
label<-as.matrix(label[[1]])

source("detec.slab2.R")
source("samr.compute.delta.table.zhang.R")
source("num_sig_zhang.R")
label[label!=1] <- 2
data=list(x=exp,y=label, geneid=rownames(exp),genenames=paste("g",as.character(1:nrow(exp)),sep=""), logged2=TRUE)
samr.obj<-samr(data, resp.type="Two class unpaired", nperms=5000)
ensemble_zhang<-num_sig_zhang(data,samr.obj,fdr=c(0.001,0.01,0.05,0.1,0.2))

write.csv(ensemble_zhang[[1]],file="result/all_sig_num.csv")
write.csv(ensemble_zhang[[2]],file="result/all_sig_gene.csv")
write.csv(samr.obj$tt,file="result/all_d_stat.csv")
write.csv(ensemble_zhang[[3]],file="result/all_cut_inf.csv")
#write.csv(samr.obj$foldchange,file="all_foldchange.csv")
}
