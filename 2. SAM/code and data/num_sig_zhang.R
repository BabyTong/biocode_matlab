num_sig_zhang<-function(data,sam.out,fdr=seq(0.01,0.1,by=0.01), FC=0)
{
	gene_num<-matrix(NA,1,length(fdr))
	cut_inf<-matrix(0,length(fdr),8)
	gene_sig<-matrix(0,dim(data$x)[[1]],length(fdr))
	for (i in 1:length(fdr)){
#browser()
		s_rst<-findFDR_samr2(data, sam.out,fdr[i], FC=FC)
#		browser()
		if (s_rst$err){
		next}
		else
		{up_index<-as.numeric(s_rst$siggenes.table$genes.up[,"Row"])-1
		lo_index<-as.numeric(s_rst$siggenes.table$genes.lo[,"Row"])-1
		gene_sig[up_index,i]<-1
		gene_sig[lo_index,i]<--1
#
		cut_inf[i,]<-s_rst$tmp[s_rst$temp,]
		gene_num[i]<-sum(abs(gene_sig[,i]))
		#gene_num[i]<-delta[1,"# called"]
		}
		}
		gene_inf<-sum_col(t(gene_sig))
		colnames(cut_inf)<-colnames(s_rst$tmp)
		return(list(gene_num=gene_num,gene_inf=gene_inf,cut_inf=cut_inf))
		}

findFDR_samr2<-function(data,sam.out,fdr=0.1,delta=NULL,FC=0,prec=6,nvals=50,...){
	require(samr)
#browser()
# added by zjf 29/5/2008
	tmp<-samr.compute.delta.table.zhang(sam.out, min=FC, dels=delta,nvals=50,...)
#browser()
#	write.csv(tmp,"D:\\tmp_zhang.csv")
#browser()
	vec.fdr<-tmp[,"90th perc FDR"]
	delta<-tmp[,"delta"]
	if(all(vec.fdr>fdr,na.rm=TRUE)){
		cat("All FDRs are larger than",fdr,"\n")
		save(list=ls(),file="error_report.Rdata")
		return(list(err=TRUE))		
	}
	if(all(vec.fdr<fdr,na.rm=TRUE)){
		cat("All FDRs are smaller than",fdr,"\n")
#		browser()
		cat("delta=",tmp[1,"delta"],"\n")
		cat("FDR=",tmp[1,"delta"],"\n")
		save(list=ls(),file="error_report.Rdata")
		return(list(err=TRUE))
	}
#	tmp2<-min(which(vec.fdr<fdr))
#added by zjf 2009-4-10£ºprocessing the situation with more than one segments lower than fdr cutoff
	tmp2 <- max(which(vec.fdr>fdr))+1
#added by zjf 2009-4-10: choosing the segments with the maximum delta region.
	delta.new<-round(seq(tmp[tmp2-1,"delta"],tmp[tmp2,"delta"],le=10),prec)
	if(all(delta==delta.new)){
		cat("It seems that the threshold is here.","\n") 
#		browser()
# added by zjf in 12/8/2007
		temp<-min(which(vec.fdr<fdr))
#		threshold<-delta.new[[1]]
		threshold<-delta.new[[temp]]
		siggenes.table<-samr.compute.siggenes.table(sam.out,threshold, data, tmp)
		return(list(err=FALSE,threshold=threshold,siggenes.table=siggenes.table,tmp=tmp,temp=temp))
               
	}
	else
		delta<-delta.new
	cat("\n","Now searching between:","\n","delta:",tmp[tmp2-1,"delta"],
		"    ","FDR:",tmp[tmp2-1,"90th perc FDR"],"\n","delta:",tmp[tmp2,"delta"],
		"    ","FDR:",tmp[tmp2,"90th perc FDR"],"\n")
	findFDR_samr2(data,sam.out,fdr=fdr,delta=delta,FC=FC,prec=prec)
	}
	
	sum_col<-function(mat){
	s_col<-t(mat)%*%rep(1,dim(mat)[[1]])
	return(s_col)
}

