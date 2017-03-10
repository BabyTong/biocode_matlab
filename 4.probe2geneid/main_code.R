setwd("C:\\Users\\Administrator\\Desktop\\shuju\\GSE20194")
expdata<-justRMA(normalize=FALSE)
write.exprs(expdata,file="probe_exp.txt")
read.table('probe_geneid.txt',header=F,sep="\t")->probeid_geneid1
read.table('probe_exp.txt',header=TRUE,sep="\t",fill=NA)->probeid_exp1
source("RMA.R")
probeid2geneid(probeid_exp1,probeid_geneid1)->profile1
write.table(profile1,'GSE6532(97)RMA.txt',sep="\t") 
