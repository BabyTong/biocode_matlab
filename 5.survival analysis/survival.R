source("http://www.bioconductor.org/biocLite.R")
biocLite("survival","release")
setwd()
class_ind<-as.matrix(read.table("ind_4922.txt",header=F,sep="\t",stringsAsFactors=FALSE))##样本分类结果
surv_info<-read.table("surv_info_4922.txt",header=T,sep="\t",stringsAsFactors=FALSE)##生存信息
muti_info<-read.table("muti_4922.txt",header=T,sep="\t",stringsAsFactors=FALSE)##所有临床信息
##注意参数，header,stringsAsFactors
library(survival)
ind<-as.vector(class_ind)

##KM曲线
surv_obj<-survfit(Surv(months, relapse) ~ ind,data=surv_info)
plot(surv_obj)

##单变量cox
summary(coxph(Surv(months, relapse) ~ ind,data=surv_info))

##多变量cox
summary(coxph(Surv(months, relapse) ~ ind+age+grade+size,data=muti_info))

##画图1
pdf("Figure.pdf",height=5,width=5)
surv_obj<-survfit(Surv(months, relapse) ~ ind,data=surv_info)
plot(surv_obj,main="GSE4922",col=c(rgb(27,74,158,max=255),rgb(238,121,66,max=255)),xlab="Recurrence-Free Survival Time(months)",ylab="Proportion Free of Recurrence",lwd=2,font.lab=2,font.axis=2,xaxt='n')
axis(side=1,at=c(0,25,50,75,100,125,150,175),labels=c(0,25,50,75,100,125,150,175),font=2)  
legend("bottomleft", legend =c("Low risk", "High risk"),col=c(rgb(27,74,158,max=255),rgb(238,121,66,max=255)),text.font=2,lty=c(1,1),lwd=2,bty="n")
legend("bottomright", legend =c("P=4.49E-03","HR=2.61(95%CI,1.31-5.19)"),text.font=2,bty="n")
box()
dev.off()

##画图3
pdf("Figure1.pdf",height=5,width=15)
layout(matrix(c(1,2,3),1,3,byrow=TRUE),c(1,1),c(1,1),TRUE)
surv_obj<-survfit(Surv(months, relapse) ~ ind,data=surv_info)
plot(surv_obj,main="GSE4922",col=c(rgb(27,74,158,max=255),rgb(238,121,66,max=255)),xlab="Recurrence-Free Survival Time(months)",ylab="Proportion Free of Recurrence",lwd=2,font.lab=2,font.axis=2,xaxt='n')
axis(side=1,at=c(0,25,50,75,100,125,150,175),labels=c(0,25,50,75,100,125,150,175),font=2)  
legend("bottomleft", legend =c("Low risk", "High risk"),col=c(rgb(27,74,158,max=255),rgb(238,121,66,max=255)),text.font=2,lty=c(1,1),lwd=2,bty="n")
legend("bottomright", legend =c("P=4.49E-03","HR=2.61(95%CI,1.31-5.19)"),text.font=2,bty="n")
surv_obj<-survfit(Surv(months, relapse) ~ ind,data=surv_info)
plot(surv_obj,main="GSE4922",col=c(rgb(27,74,158,max=255),rgb(238,121,66,max=255)),xlab="Recurrence-Free Survival Time(months)",ylab="Proportion Free of Recurrence",lwd=2,font.lab=2,font.axis=2,xaxt='n')
axis(side=1,at=c(0,25,50,75,100,125,150,175),labels=c(0,25,50,75,100,125,150,175),font=2)  
legend("bottomleft", legend =c("Low risk", "High risk"),col=c(rgb(27,74,158,max=255),rgb(238,121,66,max=255)),text.font=2,lty=c(1,1),lwd=2,bty="n")
legend("bottomright", legend =c("P=4.49E-03","HR=2.61(95%CI,1.31-5.19)"),text.font=2,bty="n")
surv_obj<-survfit(Surv(months, relapse) ~ ind,data=surv_info)
plot(surv_obj,main="GSE4922",col=c(rgb(27,74,158,max=255),rgb(238,121,66,max=255)),xlab="Recurrence-Free Survival Time(months)",ylab="Proportion Free of Recurrence",lwd=2,font.lab=2,font.axis=2,xaxt='n')
axis(side=1,at=c(0,25,50,75,100,125,150,175),labels=c(0,25,50,75,100,125,150,175),font=2)  
legend("bottomleft", legend =c("Low risk", "High risk"),col=c(rgb(27,74,158,max=255),rgb(238,121,66,max=255)),text.font=2,lty=c(1,1),lwd=2,bty="n")
legend("bottomright", legend =c("P=4.49E-03","HR=2.61(95%CI,1.31-5.19)"),text.font=2,bty="n")
box()
dev.off()



plot(surv_obj,main="GSE4922",col=c("#FFE4C4","#00C5CD"),xlab="Recurrence-Free Survival Time(months)",ylab="Proportion Free of Recurrence",lwd=2,font.lab=2,font.axis=2,xaxt='n')
