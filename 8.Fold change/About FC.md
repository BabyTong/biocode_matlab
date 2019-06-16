​	倍数法（Fold change） 是最早应用于基因表达谱数据分析的方法，用
于量化一个基因的表达变化程度：



![img](https://upload-images.jianshu.io/upload_images/1203538-e3f69b83686ad1e4.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/350/format/webp)



​	其中 x 1( i ) 和 x 2(i ) 分别表示基因 i 在两类样本中原始表达丰度的均值。 该方法通
过两类样本的基因原始表达丰度值的比较，计算出每个基因在两类样本间的
倍数差值，如果倍数差值大于预设阈值，判定为差异表达基因。虽然该方法
适用于小样本数据分析，但是没有考虑基因差异表达的统计显著性而且差异
阈值的设定相对人为，常根据研究者的经验及需要进行设定。

