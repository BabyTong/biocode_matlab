probeid2geneid=function(probeid_exp1,probeid_geneid1){
as.matrix(probeid_exp1)->probeid_exp
as.matrix(probeid_geneid1)->probeid_geneid
pid.exp=probeid_exp[,1]
index=match(probeid_geneid[,1],pid.exp)
probeid_exp=probeid_exp[index,]
raw.geneid=as.numeric(as.matrix(probeid_geneid[,2]))
drop.index=which(is.na(raw.geneid))
if(length(drop.index)>0)
{
raw.exp=probeid_exp[,-1]
mode(raw.exp)="numeric"
geneid=raw.geneid[-drop.index]
exp.matrix=raw.exp[-drop.index,]
}
else
{
raw.exp=probeid_exp[,-1];
mode(raw.exp)="numeric";
geneid=raw.geneid;
exp.matrix=raw.exp;
}
geneidfactor=factor(geneid)
exp1.matrix=apply(exp.matrix,2,function(x) tapply(x,geneidfactor,mean))

}
