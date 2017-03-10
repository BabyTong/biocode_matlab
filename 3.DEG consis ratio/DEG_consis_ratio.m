function [com_up,com_down,com_up_down,inter,ratio,p1,p2]=DEG_consis_ratio(up1,up2,down1,down2,DEG1,DEG2,Gid1,Gid2)
%[com_up,com_down,com_up_down,inter,ratio,p1,p2]=DEG_consis_ratio(up1,up2,down1,down2,DEG1,DEG2,Gid1,Gid2)
%[IN] up1 up2 :up_regulation DEG(different express Gene)
%[OUT] com_up: two group DEG commmen up_regulation
com_up=intersect(up1,up2);
com_down=intersect(down1,down2);
com_up_down=[com_up;com_down];
inter=intersect(DEG1,DEG2);
ratio=length(com_up_down)./length(inter);
Gid=intersect(Gid1,Gid2);
p1=1-binocdf(length(com_up_down)-1,length(inter),0.5);
p2=1-hygecdf(length(inter)-1,length(Gid),length(DEG1),length(DEG2));
end
