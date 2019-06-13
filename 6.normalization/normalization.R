################################
#标准化算法程序代码
################################

#include:RMA MAS5.0 dChip LVS Quantile Lowess

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#1.RMA
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#1.1主程序
justRMA<-function (..., filenames = character(0), widget = getOption("BioC")$affy$use.widgets, 
          compress = getOption("BioC")$affy$compress.cel, celfile.path = getwd(), 
          sampleNames = NULL, phenoData = NULL, description = NULL, 
          notes = "", rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE, 
          hdf5 = FALSE, hdf5FilePath = NULL, verbose = FALSE, normalize = TRUE, 
          background = TRUE, bgversion = 2, destructive = FALSE, cdfname = NULL) 
{
  l <- AllButCelsForReadAffy(..., filenames = filenames, widget = widget, 
                             celfile.path = celfile.path, sampleNames = sampleNames, 
                             phenoData = phenoData, description = description)
  ret <- just.rma(filenames = l$filenames, phenoData = l$phenoData, 
                  description = l$description, notes = notes, compress = compress, 
                  rm.mask = rm.mask, rm.outliers = rm.outliers, rm.extra = rm.extra, 
                  verbose = verbose, normalize = normalize, background = background, 
                  bgversion = bgversion, destructive = destructive, cdfname = cdfname)
  sampleNames(ret) <- l$sampleNames
  return(ret)
}

#1.2子程序
AllButCelsForReadAffy<-function (..., filenames = character(0), widget = getOption("BioC")$affy$use.widgets, 
          celfile.path = NULL, sampleNames = NULL, phenoData = NULL, 
          description = NULL) 
{
  auxnames <- unlist(as.list(substitute(list(...)))[-1])
  if (widget) {
    require(tkWidgets)
    widgetfiles <- fileBrowser(textToShow = "Choose CEL files", 
                               testFun = hasSuffix("[cC][eE][lL]"))
  }
  else {
    widgetfiles <- character(0)
  }
  if (!is.null(celfile.path)) {
    auxnames <- file.path(celfile.path, auxnames)
    filenames <- file.path(celfile.path, filenames)
  }
  filenames <- c(filenames, auxnames, widgetfiles)
  if (length(filenames) == 0) {
    if (is.null(celfile.path)) 
      celfile.path <- getwd()
    filenames <- list.celfiles(celfile.path, full.names = TRUE)
  }
  if (length(filenames) == 0) 
    stop("No cel filennames specified and no cel files in specified directory:", 
         celfile.path, "\n")
  if (is.null(sampleNames)) {
    sampleNames <- sub("^/?([^/]*/)*", "", filenames)
  }
  else {
    if (length(sampleNames) != length(filenames)) {
      warning("sampleNames not same length as filenames. Using filenames as sampleNames instead\n")
      sampleNames <- sub("^/?([^/]*/)*", "", filenames)
    }
  }
  if (is.character(phenoData)) {
    if (length(phenoData) != 1) 
      stop(sprintf("'phenoData' must be of length 1, but is %d.", 
                   length(phenoData)))
    phenoData <- read.AnnotatedDataFrame(filename = phenoData)
    sampleNames <- sampleNames(phenoData)
  }
  else if (is.data.frame(phenoData)) {
    phenoData <- as(phenoData, "AnnotatedDataFrame")
  }
  else if (is.null(phenoData)) {
    phenoData <- new("AnnotatedDataFrame", data = data.frame(sample = seq_along(sampleNames), 
                                                             row.names = sampleNames), varMetadata = data.frame(labelDescription = "arbitrary numbering", 
                                                                                                                row.names = names(pData)))
  }
  else if (!is(phenoData, "AnnotatedDataFrame")) {
    stop(sprintf("'phenoData' must be of class 'AnnotatedDataFrame', but is %s.", 
                 class(phenoData)))
  }
  if (is.character(description)) {
    description <- read.MIAME(filename = description, widget = FALSE)
  }
  else {
    if (!is(description, "MIAME")) {
      description <- new("MIAME")
    }
  }
  description@preprocessing$filenames <- filenames
  description@preprocessing$affyversion <- library(help = affy)$info[[2]][[2]][2]
  return(list(filenames = filenames, phenoData = phenoData, 
              sampleNames = sampleNames, description = description))
}

just.rma<-function (..., filenames = character(0), phenoData = new("AnnotatedDataFrame"), 
          description = NULL, notes = "", compress = getOption("BioC")$affy$compress.cel, 
          rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE, verbose = FALSE, 
          background = TRUE, normalize = TRUE, bgversion = 2, destructive = FALSE, 
          cdfname = NULL) 
{
  auxnames <- unlist(list(...))
  filenames <- c(filenames, auxnames)
  checkValidFilenames(filenames)
  n <- length(filenames)
  pdata <- pData(phenoData)
  if (dim(pdata)[1] != n) {
    warning("Incompatible phenoData object. Created a new one.\n")
    samplenames <- gsub("^/?([^/]*/)*", "", unlist(filenames))
    pdata <- data.frame(sample = 1:n, row.names = samplenames)
    phenoData <- new("AnnotatedDataFrame", data = pdata, 
                     varMetadata = data.frame(labelDescription = "arbitrary numbering", 
                                              row.names = "sample"))
  }
  else samplenames <- rownames(pdata)
  if (is.null(description)) {
    description <- new("MIAME")
    description@preprocessing$filenames <- filenames
    description@preprocessing$affyversion <- library(help = affy)$info[[2]][[2]][2]
  }
  headdetails <- read.celfile.header(filenames[[1]])
  if (is.null(cdfname)) 
    cdfname <- headdetails[[1]]
  scandates <- sapply(seq_len(length(filenames)), function(i) {
    sdate <- read.celfile.header(filenames[i], info = "full")[["ScanDate"]]
    if (is.null(sdate) || length(sdate) == 0) 
      NA_character_
    else sdate
  })
  protocol <- new("AnnotatedDataFrame", data = data.frame(ScanDate = scandates, 
                                                          row.names = sampleNames(phenoData), stringsAsFactors = FALSE), 
                  dimLabels = c("sampleNames", "sampleColumns"))
  tmp <- new("AffyBatch", cdfName = cdfname, annotation = cleancdfname(cdfname, 
                                                                       addcdf = FALSE))
  pmIndex <- pmindex(tmp)
  probenames <- rep(names(pmIndex), unlist(lapply(pmIndex, 
                                                  length)))
  pNList <- split(0:(length(probenames) - 1), probenames)
  probeintensities <- read.probematrix(filenames = filenames, 
                                       cdfname = cdfname)
  ngenes <- length(geneNames(tmp))
  exprs <- .Call("rma_c_complete", probeintensities$pm, pNList, 
                 ngenes, normalize, background, bgversion, verbose, PACKAGE = "affy")
  colnames(exprs) <- samplenames
  se.exprs <- array(NA, dim(exprs), dimnames = list(rownames(exprs), 
                                                    colnames(exprs)))
  annotation <- annotation(tmp)
  notes(description) <- notes
  new("ExpressionSet", phenoData = phenoData, protocolData = protocol, 
      annotation = annotation, experimentData = description, 
      exprs = exprs, se.exprs = se.exprs)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#2.MAS5.0
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#2.1主程序
mas5<-function (object, normalize = TRUE, sc = 500, analysis = "absolute", 
          ...) 
{
  res <- expresso(object, bgcorrect.method = "mas", pmcorrect.method = "mas", 
                  normalize = FALSE, summary.method = "mas", ...)
  if (normalize) 
    res <- affy.scalevalue.exprSet(res, sc = sc, analysis = analysis)
  return(res)
}

#2.2子程序
expresso<-function (afbatch, bg.correct = TRUE, bgcorrect.method = NULL, 
          bgcorrect.param = list(), normalize = TRUE, normalize.method = NULL, 
          normalize.param = list(), pmcorrect.method = NULL, pmcorrect.param = list(), 
          summary.method = NULL, summary.param = list(), summary.subset = NULL, 
          verbose = TRUE, widget = FALSE) 
{
  setCorrections <- function() {
    bioc.opt <- getOption("BioC")
    if (bg.correct) {
      if (is.null(bgcorrect.method)) {
        BGMethods <- bgcorrect.methods()
      }
      else {
        BGMethods <- bgcorrect.method
      }
    }
    else {
      BGMethods <- "None"
    }
    if (normalize) {
      if (is.null(normalize.method)) {
        normMethods <- normalize.methods(afbatch)
      }
      else {
        normMethods <- normalize.method
      }
    }
    else {
      normMethods <- "None"
    }
    if (is.null(pmcorrect.method)) {
      PMMethods <- pmcorrect.methods()
    }
    else {
      PMMethods <- pmcorrect.method
    }
    if (is.null(summary.method)) {
      expMethods <- generateExprSet.methods()
    }
    else {
      expMethods <- summary.method
    }
    corrections <- expressoWidget(BGMethods, normMethods, 
                                  PMMethods, expMethods, bioc.opt$affy$bgcorrect.method, 
                                  bioc.opt$affy$normalize.method, bioc.opt$affy$pmcorrect.method, 
                                  bioc.opt$affy$summary.method)
    if (!is.null(corrections)) {
      if (corrections[["BG"]] != "None") {
        bgcorrect.method <<- corrections[["BG"]]
      }
      if (corrections[["NORM"]] != "None") {
        normalize.method <<- corrections[["NORM"]]
      }
      if (corrections[["PM"]] != "None") {
        pmcorrect.method <<- corrections[["PM"]]
      }
      if (corrections[["EXP"]] != "None") {
        summary.method <<- corrections[["EXP"]]
      }
    }
    else {
      stop("Aborted by user")
    }
  }
  if (widget) {
    require(tkWidgets) || stop("library tkWidgets could not be found !")
  }
  nchips <- length(afbatch)
  if (widget) {
    setCorrections()
  }
  if (verbose) {
    if (bg.correct) {
      cat("background correction:", bgcorrect.method, "\n")
    }
    if (normalize) {
      cat("normalization:", normalize.method, "\n")
    }
    cat("PM/MM correction :", pmcorrect.method, "\n")
    cat("expression values:", summary.method, "\n")
  }
  if (bg.correct) {
    if (verbose) 
      cat("background correcting...")
    afbatch <- do.call(affy:::bg.correct, c(alist(afbatch, 
                                                  method = bgcorrect.method), bgcorrect.param))
    if (verbose) 
      cat("done.\n")
  }
  if (normalize) {
    if (verbose) 
      cat("normalizing...")
    afbatch <- do.call(affy:::normalize, c(alist(afbatch, 
                                                 normalize.method), normalize.param))
    if (verbose) 
      cat("done.\n")
  }
  eset <- computeExprSet(afbatch, summary.method = summary.method, 
                         pmcorrect.method = pmcorrect.method, ids = summary.subset, 
                         summary.param = summary.param, pmcorrect.param = pmcorrect.param)
  return(eset)
}

affy.scalevalue.exprSet<-function (eset, sc = 500, analysis = "absolute") 
{
  analysis <- match(analysis, c("absolute", "comparison"))
  if (analysis == 1) 
    nf <- 1
  else stop("sorry! comparison not implemented.")
  for (i in 1:ncol(exprs(eset))) {
    slg <- exprs(eset)[, i]
    sf <- sc/mean(slg, trim = 0.02)
    reported.value <- nf * sf * slg
    exprs(eset)[, i] <- reported.value
  }
  return(eset)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#3.dChip
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#3.1主程序

fit.li.wong<-function (data.matrix, remove.outliers = TRUE, normal.array.quantile = 0.5, 
          normal.resid.quantile = 0.9, large.threshold = 3, large.variation = 0.8, 
          outlier.fraction = 0.14, delta = 1e-06, maxit = 50, outer.maxit = 50, 
          verbose = FALSE, ...) 
{
  if (missing(data.matrix)) 
    stop("Argument data.matrix missing, with no default")
  II <- dim(data.matrix)[1]
  J <- dim(data.matrix)[2]
  if (J == 1) {
    warning("Li and Wong's algorithm is not suitable when only one probe pair")
    return(list(theta = as.vector(data.matrix), phi = 1, 
                sigma.eps = NA, sigma.theta = NA, sigma.phi = NA, 
                theta.outliers = NA, phi.outliers = NA, single.outliers = NA, 
                convergence1 = NA, convergence2 = NA, iter = NA, 
                delta = NA))
  }
  cI <- II
  cJ <- J
  theta.outliers.old <- rep(FALSE, II)
  phi.outliers.old <- rep(FALSE, J)
  single.outliers.old <- matrix(FALSE, II, J)
  theta.outliers <- theta.outliers.old
  phi.outliers <- phi.outliers.old
  single.outliers <- single.outliers.old
  flag1 <- NA
  flag2 <- NA
  if (remove.outliers) {
    flag1 <- TRUE
    flag2 <- TRUE
    original.data.matrix <- data.matrix
    change.theta <- 1
    change.phi <- 1
    change.single <- 1
    outer.iter <- 0
    while (flag1 & flag2 & change.theta + change.phi + change.single > 
           0 & outer.iter < outer.maxit) {
      outer.iter <- outer.iter + 1
      if ((outer.iter%%3 == 0 & change.theta > 0) | (outer.iter%%3 == 
                                                     1 & change.phi > 0)) {
        phi <- colMeans(data.matrix)
        c <- sqrt(cJ/sum(phi[!phi.outliers]^2))
        phi <- c * phi
        theta <- (data.matrix[, !phi.outliers, drop = FALSE] %*% 
                    phi[!phi.outliers, drop = FALSE])/cJ
        iter <- 0
        change <- 1
        theta.old <- rep(0, II)
        while (change > delta & iter < maxit) {
          iter <- iter + 1
          phi <- t(data.matrix[!theta.outliers, , drop = FALSE]) %*% 
            theta[!theta.outliers, drop = FALSE]
          c <- sqrt(cJ/sum(phi[!phi.outliers, drop = FALSE]^2))
          phi <- c * phi
          theta <- (data.matrix[, !phi.outliers, drop = FALSE] %*% 
                      phi[!phi.outliers, drop = FALSE])/cJ
          change <- max(abs(theta[!theta.outliers] - 
                              theta.old[!theta.outliers]))
          if (verbose) 
            cat(paste("Outlier iteration:", outer.iter, 
                      "estimation iteration:", iter, "chage=", 
                      change, "\n"))
          theta.old <- theta
        }
        if (iter >= maxit) {
          warning(paste("No convergence in inner loop after", 
                        iter, "in outerler tieration", outer.iter, 
                        "\n"))
          flag1 <- FALSE
        }
        if (mean(phi[!phi.outliers] < 0) > 0.5) {
          theta <- -theta
          phi <- -phi
        }
        theta <- as.vector(theta)
        phi <- as.vector(phi)
        data.matrixhat <- outer(theta, phi)
        resid <- data.matrix - data.matrixhat
      }
      if (outer.iter%%3 == 1) {
        single.outliers <- resid > large.threshold * 
          quantile(abs(resid), normal.resid.quantile)
        single.outliers[rowSums(single.outliers) > outlier.fraction * 
                          cJ, ] <- rep(FALSE, J)
        single.outliers[, colSums(single.outliers) > 
                          outlier.fraction * cI] <- rep(FALSE, II)
        data.matrix[single.outliers] <- data.matrixhat[single.outliers]
        data.matrix[!single.outliers] <- original.data.matrix[!single.outliers]
        change.single <- sum(abs(single.outliers.old - 
                                   single.outliers))
        single.outliers.old <- single.outliers
      }
      else {
        sigma.theta <- sqrt(rowSums(resid[, !phi.outliers, 
                                          drop = FALSE]^2)/(cJ - 1))
        sigma.phi <- sqrt(colSums(resid[!theta.outliers, 
                                        , drop = FALSE]^2)/(cI - 1))
        if (outer.iter%%3 == 2) {
          theta.outliers <- sigma.theta > large.threshold * 
            quantile(sigma.theta, normal.array.quantile) | 
            theta^2/sum(theta^2) > large.variation
          cI <- sum(!theta.outliers)
          if (cI < 3) {
            warning("No convergence achieved, too many outliers")
            flag2 <- FALSE
          }
          single.outliers[theta.outliers, ] <- rep(FALSE, 
                                                   J)
          data.matrix[single.outliers] <- data.matrixhat[single.outliers]
          data.matrix[!single.outliers] <- original.data.matrix[!single.outliers]
          change.theta <- sum(abs(theta.outliers.old - 
                                    theta.outliers))
          change.single <- sum(abs(single.outliers.old - 
                                     single.outliers))
          theta.outliers.old <- theta.outliers
        }
        else {
          phi.outliers <- sigma.phi > large.threshold * 
            quantile(sigma.phi, normal.array.quantile) | 
            phi^2/sum(phi^2) > large.variation | phi < 
            0
          cJ <- sum(!phi.outliers)
          if (cJ < 3) {
            warning("No convergence achieved, too many outliers")
            flag2 <- FALSE
          }
          single.outliers[, phi.outliers] <- rep(FALSE, 
                                                 II)
          data.matrix[single.outliers] <- data.matrixhat[single.outliers]
          data.matrix[!single.outliers] <- original.data.matrix[!single.outliers]
          change.phi <- sum(abs(phi.outliers.old - phi.outliers))
          change.single <- sum(abs(single.outliers.old - 
                                     single.outliers))
          phi.outliers.old <- phi.outliers
        }
      }
      if (verbose) {
        cat("chips used=", cI, ", probes used=", cJ, 
            ", single outler=", sum(single.outliers), "\n")
        cat("Number of changes: single=", change.single, 
            ", theta=", change.theta, ", phi=", change.phi, 
            "\n", sep = "")
      }
    }
    if (outer.iter >= outer.maxit) {
      warning("No convergence achieved in outlier loop\n")
      flag2 <- FALSE
    }
    all.outliers <- outer(theta.outliers, phi.outliers, FUN = "|") | 
      single.outliers
    sigma <- sqrt(sum(resid[!all.outliers]^2)/sum(!all.outliers))
    sigma.theta <- sqrt(rowSums(resid[, !phi.outliers, drop = FALSE]^2)/(cJ - 
                                                                           1))
    sigma.phi <- sqrt(colSums(resid[!theta.outliers, , drop = FALSE]^2)/(cI - 
                                                                           1))
  }
  else {
    flag1 <- TRUE
    phi <- colMeans(data.matrix)
    c <- sqrt(J/sum(phi^2))
    phi <- c * phi
    theta <- (data.matrix %*% phi)/J
    iter <- 0
    change <- 1
    theta.old <- rep(0, II)
    while (change > delta & iter < maxit) {
      iter <- iter + 1
      phi <- t(data.matrix) %*% theta
      c <- sqrt(J/sum(phi^2))
      phi <- c * phi
      theta <- (data.matrix %*% phi)/J
      change <- max(abs(theta - theta.old))
      if (verbose) 
        cat(paste("Iteration:", iter, "chage=", change, 
                  "\n"))
      theta.old <- theta
    }
    if (iter >= maxit) {
      warning(paste("No convergence after", iter, "iterations.\n"))
      flag1 <- FALSE
    }
    if (mean(phi[!phi.outliers] < 0) > 0.5) {
      theta <- -theta
      phi <- -phi
    }
    theta <- as.vector(theta)
    phi <- as.vector(phi)
    data.matrixhat <- outer(theta, phi)
    sigma.theta <- sqrt(rowSums((data.matrix - data.matrixhat)^2)/(J - 
                                                                     1))
    sigma.phi <- sqrt(colSums((data.matrix - data.matrixhat)^2)/(II - 
                                                                   1))
    sigma <- sqrt(sum((data.matrix - data.matrixhat)^2)/(II * 
                                                           J))
  }
  return(list(theta = theta, phi = phi, sigma.eps = sigma, 
              sigma.theta = sigma.theta, sigma.phi = sigma.phi, theta.outliers = theta.outliers, 
              phi.outliers = phi.outliers, single.outliers = single.outliers, 
              convergence1 = flag1, convergence2 = flag2, iter = iter, 
              delta = change))
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#4.LVS
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#4.1主程序

normalize.lvs<-function (object, ref.fun = c("median", "mean"), lvs.id, use.loess = FALSE, 
          ...) 
{
  if (missing(lvs.id)) 
    stop("Must provide ids for LVS genes")
  if (!is.logical(lvs.id) && !is.numeric(lvs.id)) 
    stop("lvs.id must be a vector either logical or numeric")
  nk <- ncol(object)
  ref.fun <- match.arg(ref.fun)
  ref.data <- switch(ref.fun, median = rowMedians(object), 
                     rowMeans(object))
  out <- NULL
  for (i in 1:nk) {
    if (use.loess) {
      sm <- loess(object[lvs.id, i] ~ ref.data[lvs.id])
      a = approx(sm$fit, sm$x, xout = object[, i], rule = 2)$y
    }
    else {
      sm = smooth.spline(y = object[lvs.id, i], x = ref.data[lvs.id])
      a = approx(sm$y, sm$x, xout = object[, i], rule = 2)$y
    }
    out = cbind(out, a)
  }
  colnames(out) <- colnames(object)
  rownames(out) <- rownames(object)
  return(out)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#5.Quantile
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#5.1主程序

normalize.AffyBatch.quantiles<-function (abatch, type = c("separate", "pmonly", "mmonly", "together")) 
{
  type <- match.arg(type)
  if ((type == "pmonly") | (type == "separate")) {
    pms <- unlist(pmindex(abatch))
    noNA <- rowSums(is.na(intensity(abatch)[pms, , drop = FALSE])) == 
      0
    pms <- pms[noNA]
    intensity(abatch)[pms, ] <- normalize.quantiles(intensity(abatch)[pms, 
                                                                      , drop = FALSE], copy = FALSE)
  }
  if ((type == "mmonly") | (type == "separate")) {
    mms <- unlist(mmindex(abatch))
    noNA <- rowSums(is.na(intensity(abatch)[mms, , drop = FALSE])) == 
      0
    mms <- mms[noNA]
    intensity(abatch)[mms, ] <- normalize.quantiles(intensity(abatch)[mms, 
                                                                      , drop = FALSE], copy = FALSE)
  }
  if (type == "together") {
    pms <- unlist(indexProbes(abatch, "both"))
    intensity(abatch)[pms, ] <- normalize.quantiles(intensity(abatch)[pms, 
                                                                      , drop = FALSE], copy = FALSE)
  }
  return(abatch)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#6.Lowess
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#6.1主程序

loess.normalize<-function (mat, subset = sample(1:(dim(mat)[2]), 5000), epsilon = 10^-2, 
          maxit = 1, log.it = TRUE, verbose = TRUE, span = 2/3, family.loess = "symmetric") 
{
  J <- dim(mat)[2]
  II <- dim(mat)[1]
  newData <- mat
  if (log.it) {
    mat <- log2(mat)
    newData <- log2(newData)
  }
  change <- epsilon + 1
  fs <- matrix(0, II, J)
  iter <- 0
  w <- c(0, rep(1, length(subset)), 0)
  while (iter < maxit) {
    iter <- iter + 1
    means <- matrix(0, II, J)
    for (j in 1:(J - 1)) {
      for (k in (j + 1):J) {
        y <- newData[, j] - newData[, k]
        x <- (newData[, j] + newData[, k])/2
        index <- c(order(x)[1], subset, order(-x)[1])
        xx <- x[index]
        yy <- y[index]
        aux <- loess(yy ~ xx, span = span, degree = 1, 
                     weights = w, family = family.loess)
        aux <- predict(aux, data.frame(xx = x))/J
        means[, j] <- means[, j] + aux
        means[, k] <- means[, k] - aux
        if (verbose) 
          cat("Done with", j, "vs", k, " in iteration ", 
              iter, "\n")
      }
    }
    fs <- fs + means
    newData <- mat - fs
    change <- max(colMeans((means[subset, ])^2))
    if (verbose) 
      cat(iter, change, "\n")
    oldfs <- fs
  }
  if (change > epsilon & maxit > 1) 
    warning(paste("No convergence after", maxit, "iterations.\n"))
  if (log.it) 
    return(2^newData)
  else return(newData)
}


