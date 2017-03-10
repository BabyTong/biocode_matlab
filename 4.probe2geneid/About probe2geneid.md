# About probe2geneid

**1.Useage**

```R
setwd("C:\\Users\\Administrator\\Desktop\\shuju\\GSE20194")
expdata<-justRMA(normalize=FALSE)
write.exprs(expdata,file="probe_exp.txt")
read.table('probe_geneid.txt',header=F,sep="\t")->probeid_geneid1
read.table('probe_exp.txt',header=TRUE,sep="\t",fill=NA)->probeid_exp1
source("RMA.R")
probeid2geneid(probeid_exp1,probeid_geneid1)->profile1
write.table(profile1,'GSE6532(97)RMA.txt',sep="\t")
```





**2.描述：**

   我们下载GEO上处理过的数据时，他通常是 probe_exp 矩阵。这时候，我们需要把他转换成GeneId_exp矩阵。

注：probe_geneis.txt指的是： 平台文件 （只保留 probe &gene id 两列）

注：一个探针对应多个基因则删除这里的多个基因，多个探针对应同一个基因则保留该基因的表达均值

​			