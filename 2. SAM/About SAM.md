# About SAM

1.useage:

```R
#sam_demo.R
setwd("F:\\程序\\samr\\")
source("http://bioconductor.org/biocLite.R")
biocLite("samr")
setwd("E:\\大创\\程序\\SAM")
library(samr)
source("integrate.r")
integrate()
```

```matlab
%sam.m
sam_deg=[geneid,allsiggene];
sam_deg(:,3)=(sam_deg(:,2)>0);
sam_deg(:,2)=abs(sam_deg(:,2));
sam001=sam_deg(sam_deg(:,2)==5|sam_deg(:,2)==4,:);
sam_F3_005=sam_deg(sam_deg(:,2)==5|sam_deg(:,2)==4|sam_deg(:,2)==3,:);
sam_F3_01=sam_deg(sam_deg(:,2)==5|sam_deg(:,2)==4|sam_deg(:,2)==3|sam_deg(:,2)==2,:);
sam_deg_02=sam_deg(sam_deg(:,2)==5|sam_deg(:,2)==4|sam_deg(:,2)==3|sam_deg(:,2)==2|sam_deg(:,2)==1,:);
```

*注：code is write by R* 



2.输入输出说明：

代码目录底下创建一个叫做input的文件夹，包含：

exp.csv geneid.csv label.csv 

sam.m输入输出说明：略



**3.描述：**

2001年Tusher等提出了SAM（Significance analysis of microarrays）算法。其计算公式如下：

![img](file:///C:\Users\ADMINI~1\AppData\Local\Temp\ksohtml\wps8E5D.tmp.png)

​        其中,![img](file:///C:\Users\ADMINI~1\AppData\Local\Temp\ksohtml\wpsB98.tmp.png)为样本残差标准误的校正值。在t统计量的分母中增加了一个较小的正值 ，减小了t检验的不稳定性，避免了使在两类样本中表达水平变化较小的基因因为具有较小的标准误而被误判为差异表达基因。SAM采用permutation算法估计错误发现率(false discovery rate, FDR),从而控制多重检验错误率，降低了结果的假阳性率。但由于SAM算法是以t检验为基础，但它依旧存在与t检验相似的问题：偏向于识别在两类条件下表达水平低但倍数变化大的差异基因。

