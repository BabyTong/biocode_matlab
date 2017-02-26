function [Sig_possition_all,Sig_possition_up,Sig_possition_down] = ttestfuntion(group1,group2,Gid,FDR)
%[IN] group1 group2 is your different kind of sample. 
%[IN]FDR is Significant Level that you cut
%[out] Sig_possition_all is all DGE(different gene express) . 
%[out] Sig_possition_up 差异上调基因  Sig_possition_down 差异下调基因
beta_N=group1;
beta_T=group2;
[~,p,~,tstates]=ttest2(beta_T',beta_N');
q = mafdr(p,'BHFDR',true);
Sig_possition_all = Gid(q <= FDR,1);
Sig_possition_up = Gid(q <= FDR & tstates.tstat>0,1);
Sig_possition_down = Gid(q <= FDR & tstates.tstat<0,1);
end

