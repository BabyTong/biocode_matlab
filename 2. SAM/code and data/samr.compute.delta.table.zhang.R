samr.compute.delta.table.zhang<-function (samr.obj, min.foldchange = 0, dels = NULL, nvals = 50) 
{
    lmax = sqrt(max(abs(sort(samr.obj$tt) - samr.obj$evo)))
    if (is.null(dels)) {
        dels = (seq(0, lmax, length = nvals)^2)
    }

    col = matrix(1, nrow = length(samr.obj$evo), ncol = nvals)
    ttstar0 <- samr.obj$ttstar0
    tt <- samr.obj$tt
    n <- samr.obj$n
    evo <- samr.obj$evo
    nsim <- ncol(ttstar0)
    res1 <- NULL
    ttstar0copy<-matrix(NA,nrow(samr.obj$ttstar0),ncol(samr.obj$ttstar0))
    foldchange.cond.up = matrix(T, nrow = nrow(samr.obj$ttstar), 
        ncol = ncol(samr.obj$ttstar))
    foldchange.cond.lo = matrix(T, nrow = nrow(samr.obj$ttstar), 
        ncol = ncol(samr.obj$ttstar))
    if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
        0)) {
        foldchange.cond.up = samr.obj$foldchange.star >= min.foldchange
        foldchange.cond.lo = samr.obj$foldchange.star <= 1/min.foldchange
    }
    cutup = rep(NA, length(dels))
    cutlow = rep(NA, length(dels))
    g2 = rep(NA, length(dels))
    errup = matrix(NA, ncol = length(dels), nrow = ncol(samr.obj$ttstar0))
    errlow = matrix(NA, ncol = length(dels), nrow = ncol(samr.obj$ttstar0))
    for (ii in 1:length(dels)) {
        cat(ii, fill = TRUE)
        ttt <- detec.slab2(samr.obj, dels[ii], min.foldchange)
#browser()
        cutup[ii] <- 1e+10
        if (length(ttt$pup > 0)) {
            cutup[ii] <- min(samr.obj$tt[ttt$pup])
        }
        cutlow[ii] <- -1e+10
        if (length(ttt$plow) > 0) {
            cutlow[ii] <- max(samr.obj$tt[ttt$plow])
        }
        g2[ii] = sumlengths2(ttt)
#browser()
# added by zjf in 29/5/2008
	 g2rows<-union(ttt$pup,ttt$plow)
	 ttstar0copy<-samr.obj$ttstar0
        ttstar0copy[g2rows,]<-0
#	 nong2rows<-setdiff(c(1:samr.obj$n),g2rows)
        errup[, ii] = colSums(ttstar0copy> cutup[ii] & 
            foldchange.cond.up)
        errlow[, ii] = colSums(ttstar0copy < cutlow[ii] & 
            foldchange.cond.lo)
    }
# added by zjf in 29/5/2008
#browser()
    s <- sqrt(apply(errup, 2, var)/nsim + apply(errlow, 2, var)/nsim)
    gmed <- apply(errup + errlow, 2, median)
    gmed_2 <- (samr.obj$n/(rep(samr.obj$n,length(g2))-g2))*gmed
    g90 = apply(errup + errlow, 2, quantile, 0.9)
    g90_2 <- (samr.obj$n/(rep(samr.obj$n,length(g2))-g2))*g90

# added by zjf in 29/5/2008
    res1 <- cbind(samr.obj$pi0*gmed_2, samr.obj$pi0*g90_2, g2, samr.obj$pi0*gmed_2/g2,samr.obj$pi0*g90_2/g2, cutlow,cutup)
# added by zjf in 29/5/2008

    res1 <- cbind(dels, res1)
    dimnames(res1) <- list(NULL, c("delta", "# med false pos", 
        "90th perc false pos", "# called", "median FDR", "90th perc FDR", 
        "cutlo", "cuthi"))
    return(res1)
}

sumlengths2<-function (aa) 
{
    length(aa$pl) + length(aa$pu)
}