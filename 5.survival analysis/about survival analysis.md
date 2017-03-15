# about survival analysis

**1.useage**

```R
# see survival.R
```

**2.描述**

输入：三个文件

**surv_info_4922.txt 包含：**months随访时间和relapse生存结局 两列

**muti_4922.txt包含：**months relapse age grade size  五列
**ind.txt ：**是自己所创建的分类器分出来的结果

除了months 这个指标，其他按照设定的规则全部0 1化

> months随访时间
>
> relapse生存结局
>
> age	 >55 <55
>
> grade 阶段
>
> size  
>
> ind 是自己所创建的分类器分出来的结果

[生存分析代码例子](http://www.bio-info-club.com/?p=246)链接里有输入的解释和输出的展示。

