foldChange<-function(inData,classLabel)
{
  #Calculating all probes' FC value
  sampleIdsCase<-which(classLabel==0);#0 tumer
  sampleIdsControl<-which(classLabel==1);#1 normal
  probeFC<-rep(0,nrow(inData))
  for(i in 1:nrow(inData))
  {
    probeFC[i]<-mean(as.numeric(as.vector(inData[i,sampleIdsCase])))/mean(as.numeric(as.vector(inData[i,sampleIdsControl])));
  }
  probeFC<-log(probeFC,base=2);
  result<-probeFC;
  return(result)
}