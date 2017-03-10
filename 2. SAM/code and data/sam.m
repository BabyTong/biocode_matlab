sam_deg=[geneid,allsiggene];
sam_deg(:,3)=(sam_deg(:,2)>0);
sam_deg(:,2)=abs(sam_deg(:,2));
sam001=sam_deg(sam_deg(:,2)==5|sam_deg(:,2)==4,:);
sam_F3_005=sam_deg(sam_deg(:,2)==5|sam_deg(:,2)==4|sam_deg(:,2)==3,:);
sam_F3_01=sam_deg(sam_deg(:,2)==5|sam_deg(:,2)==4|sam_deg(:,2)==3|sam_deg(:,2)==2,:);
sam_deg_02=sam_deg(sam_deg(:,2)==5|sam_deg(:,2)==4|sam_deg(:,2)==3|sam_deg(:,2)==2|sam_deg(:,2)==1,:);
