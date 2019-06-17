
########################################
# t-test R
########################################

#sub function

calT <- function(inData, classLabel) {
  m <- length(ind1 <- which(classLabel == 1))
  n <- length(ind2 <- which(classLabel == 0))
  inData1 <- inData[, ind1,drop=FALSE]
  inData2 <- inData[, ind2,drop=FALSE]
  rmean1 <- rowMeans(inData1)
  rmean2 <- rowMeans(inData2)
  ss1 <- rowSums((inData1 - rmean1)^2)
  ss2 <- rowSums((inData2 - rmean2)^2)
  tt <- (m + n - 2)^0.5 * (rmean2 - rmean1)/((1/m + 1/n) *                                             (ss1 + ss2))^0.5
  return(list(T = tt, df = m + n - 2))
}


#main function

#exp_data A data frame, the expression profile to calculate p-value for each gene, the rownames should be the symbol of genes.
#label A vector of 0/1s, indicating the class of samples in the expression profile, 0 represents case, 1 represents control.


t<-calT(exp_data,label)
p.t<-pt(abs(t$T),df = t$df,lower.tail = FALSE)
fdr.t<-p.adjust(p.t,method="fdr")