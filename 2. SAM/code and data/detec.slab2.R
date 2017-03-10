detec.slab2<-function (samr.obj, del, min.foldchange) 
{
    n <- length(samr.obj$tt)
    tt <- samr.obj$tt
    evo <- samr.obj$evo
    numer <- samr.obj$tt * (samr.obj$sd + samr.obj$s0)
    tag <- order(tt)
    pup <- NULL
    foldchange.cond.up = rep(T, length(evo))
    foldchange.cond.lo = rep(T, length(evo))
    if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
        0)) {
        foldchange.cond.up = samr.obj$foldchange >= min.foldchange
        foldchange.cond.lo = samr.obj$foldchange <= 1/min.foldchange
    }
    o1 <- (1:n)[(tt[tag] - evo > del) & evo > 0]
    if (length(o1) > 0) {
        o1 <- o1[1]
        o11 <- o1:n
        o111 <- rep(F, n)
        o111[tag][o11] <- T
        pup <- (1:n)[o111 & foldchange.cond.up]
    }
    plow <- NULL
    o2 <- (1:n)[(evo - tt[tag] > del) & evo < 0]
    if (length(o2) > 0) {
        o2 <- o2[length(o2)]
        o22 <- 1:o2
        o222 <- rep(F, n)
        o222[tag][o22] <- T
        plow <- (1:n)[o222 & foldchange.cond.lo]
    }
    return(list(plow = plow, pup = pup))
}