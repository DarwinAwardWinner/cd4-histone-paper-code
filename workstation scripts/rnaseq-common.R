source("common.R")

library(rtracklayer)
library(partitions)
library(biomaRt)
library(Biobase)
library(splines)

library(arrayQualityMetrics)
library(edgeR)
library(DESeq)
library(baySeq)
## library(DSS)
library(cqn)

read.counts <- function(file, assayname="counts") {
  within(list(), {
    sexp <- readRDS(file)
    stopifnot(is(sexp, "SummarizedExperiment"))
    ## colData might be a BamFileList with mcols
    if (is(colData(sexp)[[1]], "RsamtoolsFileList")) {
      colData(sexp) <- mcols(colData(sexp)[[1]])
    }
    expdata <- as(colData(sexp), "data.frame")
    counts <- assays(sexp)[[assayname]]
    gene.annot <-
      data.frame(row.names=rownames(sexp),
                 llply(as.list(mcols(rowData(sexp))),
                       function (col) {
                         if (is(col, "CharacterList")) {
                           clist.to.character(col)
                         } else {
                           as.vector(col)
                         }
                       }, .parallel=TRUE))
  })
}

generate.factor.groups <- function(v, include.null=FALSE, name.format="%s") {
  u <- droplevels(as.factor(unique(v)))
  stopifnot(length(u) >= 1)
  parts <- setparts(length(u))
  if (!include.null) {
    parts <- parts[,-1, drop=FALSE]
  }
  stopifnot(ncol(parts) >= 1)
  x <- list()
  for(i in 1:ncol(parts)) {
    partNums <- parts[,i]
    partSets <- split(u, partNums)
    partNames <- factor(foreach(s=partSets, .combine=c) %do% str_c(sort(as.character(s)), collapse="."))
    partAssignments <- partNames[partNums]
    names(partAssignments) <- u
    partLabel <- sprintf(name.format, str_c(sort(partNames), collapse="_vs_"))
    x[[partLabel]] <- partAssignments[v]
  }
  x
}

generate.all.groups <- function(df, include.null=TRUE, null.name="NDE") {
  ## Ensure everything is factors, and only keep factors with multiple
  ## levels.
  df <- data.frame(Filter(function(x) nlevels(x) > 1, Map(factor, df)))
  factor.groups <- do.call(c, setNames(Map(generate.factor.groups, df,
                                           name.format=sprintf("%s:%%s", names(df))),
                                       NULL))
  if (include.null) {
    null.group <- factor(rep("X", length(factor.groups[[1]])))
    factor.groups <- c(setNames(list(null.group), null.name),
                       factor.groups)
  }
  ## Eliminate redundant groupings
  numeric.groups <- Map(function(v) as.numeric(factor(v, levels=as.character(unique(v)))), factor.groups)
  factor.groups[!duplicated(numeric.groups)]
}

read.htseq.counts <- function(f) {
  x <- read.table(f, header=FALSE, sep="\t")
  names(x) <- c("ID", "count")
  ## Discard the last 5 lines
  exception.rows <- (nrow(x)-4):nrow(x)
  attr(x, "exception_counts") <- x[exception.rows,]
  x <- x[-exception.rows,]
  x
}

setAs(from="CountDataSet", to="ExpressionSet", def=function(from) {
    if (all(is.na(sizeFactors(from)))) {
        from <- estimateSizeFactors(from)
    }
    from <- estimateDispersions(from, method="blind", fitType="parametric")
    DESeq::varianceStabilizingTransformation(from)
})

plotMDS.CountDataSet <- function(x, ...) plotMDS(as(x, "ExpressionSet"), ...)

## The following two functions allow the use to edgeR's infrastructure
## to execute the DESeq method. Instead of glmFit(dge, design), use
## glmFit.CountDataSet(cds, design), then continue as with a normal
## edgeR analysis.
getOffset.CountDataSet <- function(y) {
  if (any(is.na(sizeFactors(y))))
    stop("Call estimateSizeFactors first")
  log(sizeFactors(y)) - mean(log(sizeFactors(y))) + mean(log(colSums(counts(y))))
}

glmFit.CountDataSet <- function (y, design = NULL, dispersion = NULL, offset = NULL,
                                 weights = NULL, lib.size = NULL, start = NULL, method = "auto",
                                 ...)
{
  stopifnot(is(y, "CountDataSet"))
  if (is.null(dispersion)) {
    if ("disp_pooled" %in% colnames(fData(y)))
      dispersion <- fData(y)$disp_pooled
    else if ("disp_blind" %in% colnames(fData(y))) {
      if (fitInfo(y, "blind")$sharingMode != "fit-only")
        warning("You have used 'method=\"blind\"' in estimateDispersion without also setting 'sharingMode=\"fit-only\"'. This will not yield useful results.")
      dispersion <- fData(y)$disp_blind
    }
    else stop("Call 'estimateDispersions' with 'method=\"pooled\"' (or 'blind') first.")
  }
  if (is.null(offset) && is.null(lib.size))
    offset <- getOffset.CountDataSet(y)
  ## UGLY HACK: DESeq can produce infinite dispersion estimates, which
  ## cause errors in glmFit. To fix this, we replace infinite
  ## dispersion estimates with the maximum representable floating
  ## point value, which should always result in a PValue of 1.
  infdisp <- !is.finite(dispersion)
  dispersion[infdisp] <- .Machine$double.xmax
  fit <- glmFit.default(y = counts(y), design = design, dispersion = dispersion,
                        offset = offset, weights = weights, lib.size = lib.size,
                        start = start, method = method, ...)
  ## Now set the dispersions back to infinite in the resulting fit object.
  fit$dispersion[infdisp] <- +Inf
  dimnames(fit$coefficients) <- list(rownames(counts(y)), colnames(design))
  fit$counts <- counts(y)
  fit$samples <- pData(y)[-1]
  fit$genes <- fData(y)[setdiff(names(fData(y)), c("disp_blind", "disp_pooled"))]
  new("DGEGLM", fit)
}

calcMA.CD <- function (cD, samplesA, samplesB, normaliseData = TRUE, log.base=2)
{
    if (is.character(samplesA)) {
        Asamps <- which(as.character(cD@replicates) %in% samplesA)
        if (!all(samplesA %in% cD@replicates))
            Asamps <- c(Asamps, which(colnames(cD@data) %in%
                samplesA[!(samplesA %in% as.character(cD@replicates))]))
        if (!all(samplesA %in% c(colnames(cD@data), as.character(cD@replicates))))
            warning("Some members of 'samplesA' were not found!")
        samplesA <- Asamps
    }
    if (length(samplesA) == 0)
        stop("Can't find any data for sample set A!")
    if (is.character(samplesB)) {
        Bsamps <- which(as.character(cD@replicates) %in% samplesB)
        if (!all(samplesB %in% cD@replicates))
            Bsamps <- c(Bsamps, which(colnames(cD@data) %in%
                samplesB[!(samplesB %in% as.character(cD@replicates))]))
        if (!all(samplesB %in% c(colnames(cD@data), as.character(cD@replicates))))
            warning("Some members of 'samplesB' were not found!")
        samplesB <- Bsamps
    }
    if (length(samplesB) == 0)
        stop("Can't find any data for sample set B!")
    if (!inherits(cD, what = "countData"))
        stop("variable 'cD' must be of or descend from class 'countData'")
    Adata <- cD@data[, samplesA]
    Bdata <- cD@data[, samplesB]
    if (normaliseData) {
        Adata <- Adata/cD@libsizes[samplesA] * mean(cD@libsizes[c(samplesA,
            samplesB)])
        Bdata <- Bdata/cD@libsizes[samplesB] * mean(cD@libsizes[c(samplesA,
            samplesB)])
    }
    if (nrow(cD@seglens) > 0)
        if (ncol(cD@seglens) == 1) {
            Adata <- Adata/cD@seglens[, 1]
            Bdata <- Bdata/cD@seglens[, 1]
        }
        else {
            Adata <- Adata/cD@seglens[, samplesA]
            Bdata <- Bdata/cD@seglens[, samplesB]
        }
    Adata <- colSums(t(Adata))/length(samplesA)
    Bdata <- colSums(t(Bdata))/length(samplesB)
    Azeros <- which(Adata == 0)
    Bzeros <- which(Bdata == 0)
    nonzeros <- which(Adata != 0 & Bdata != 0)
    bothzeros <- which(Adata == 0 & Bdata == 0)
    infRatio <- ceiling(max(abs((log(Adata, log.base) - log(Bdata, log.base))[nonzeros]),
        na.rm = TRUE))
    M <- log(Adata, log.base) - log(Bdata, log.base)
    M[Azeros] <- -infRatio - 2
    M[Bzeros] <- infRatio + 2
    M[bothzeros] <- NA
    A <- (log(Adata, log.base) + log(Bdata, log.base))/2
    A[Azeros] <- log(Bdata[Azeros], log.base)
    A[Bzeros] <- log(Adata[Bzeros], log.base)
    A[bothzeros] <- NA
    data.frame(M=M, A=A, row.names=row.names(cD@data))
}

calcFPKM <- function(counts, libsizes, seglens) {
  stopifnot(length(libsizes) == ncol(counts))
  libsizes <- matrix(libsizes, ncol=ncol(counts), nrow=nrow(counts), byrow=TRUE)
  ## Seglens can have either a single column (which could be a vector)
  ## or one column per sample
  if (is.vector(seglens)) {
    stopifnot(length(seglens) == nrow(counts))
    seglens <- matrix(seglens, ncol=ncol(counts), nrow=nrow(counts))
  } else {
    stopifnot(nrow(seglens) == nrow(counts))
    stopifnot(ncol(seglens) %in% c(1, ncol(counts)))
    seglens <- matrix(seglens, ncol=ncol(counts), nrow=nrow(counts))
  }
  counts * 1e9 / (libsizes * seglens)
}

calcFPKM.CD <- function(cD) {
  calcFPKM(cD@data, cD@libsizes,
           if (nrow(cD@seglens) > 0) cD@seglens else rep(1000, nrow(cD@data)))
}

calcFPKM.byGroup <- function(cD, group.factor) {
  group.factor <- droplevels(as.factor(group.factor))
  counts.bygroup <- splat(cbind)(llply(levels(group.factor), function (glevel) {
    rowSums(cD@data[,group.factor == glevel, drop=FALSE])
  }))
  colnames(counts.bygroup) <- levels(group.factor)
  libsizes.bygroup <- unlist(llply(levels(group.factor), function (glevel) {
    sum(cD@libsizes[group.factor == glevel])
  }))
  names(libsizes.bygroup) <- levels(group.factor)
  calcFPKM(counts.bygroup, libsizes.bygroup,
           if (nrow(cD@seglens) > 0) cD@seglens else rep(1000, nrow(cD@data)))
}

topCounts.multigroup <- function (cD, groups, decreasing = TRUE, number = 10, likelihood,
                                  FDR, normaliseData = FALSE) {
  if (!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")
  if (nrow(cD@posteriors) == 0)
    stop("The '@posteriors' slot of cD is empty!")
  if (normaliseData)
    data <- round(t(t(cD@data)/cD@libsizes) * (prod(cD@libsizes))^(1/length(cD@libsizes)))
  else data <- cD@data
  if (is.character(groups))
    groups <- pmatch(groups, names(cD@groups))
  if (nrow(cD@annotation) == 0)
    annotation <- data.frame(rowID = paste("row", 1:nrow(cD),
                               sep = "_"))
  else annotation <- cD@annotation
  if (class(cD) == "lociData")
    annotation <- cbind(data.frame(chr = as.character(seqnames(cD@coordinates)),
                                   start = as.numeric(start(cD@coordinates)), end = as.numeric(end(cD@coordinates))),
                        annotation)
  else annotation <- annotation
  if (is.null(groups)) {
    if (length(cD@nullPosts) == 0)
      stop("The '@nullPosts' slot of cD is empty - you can't use 'groups = NULL'.")
    likes <- cD@nullPosts
  }
  else likes <- log(rowSums(exp(cD@posteriors[, groups, drop = FALSE])))
  if (!missing(likelihood))
    cutNumber <- sum(likes > log(likelihood))
  if (missing(likelihood) & !missing(FDR))
    cutNumber <- sum(cumsum(1 - exp(sort(likes, decreasing = decreasing)))/1:length(likes) <
                     FDR)
  if (!missing(likelihood) | !missing(FDR))
    if (cutNumber == 0)
      warning("No features were found using the cutoffs for likelihood or FDR specified; using the 'number' argument instead")
    else number <- cutNumber
  selTags <- order(likes, decreasing = decreasing)[1:number]
  topTags <- data.frame(annotation[selTags, , drop = FALSE],
                        data[selTags, , drop = FALSE], Likelihood = exp(likes[selTags]),
                        FDR = cumsum(1 - exp(likes[selTags]))/1:number)
  rownames(topTags) <- rownames(cD@data)[selTags]
  topTags
}

topCounts.by.group.type <- function(cD, group.factor, ...) {
  group.factor <- droplevels(as.factor(group.factor))
  setNames(nm=levels(group.factor),
           llply(levels(group.factor), function(group.type) {
             groups <- which(group.factor == group.type)
             topCounts.multigroup(cD, groups, ...)
           }, .parallel=TRUE))
}

calcFC <- function(data, input.log=FALSE, output.log=input.log, strip.prefix=NULL) {
  data <- as.matrix(data)
  if (is.null(colnames(data))) {
    colnames(data) <- make.names(as.character(1:ncol(data)))
  } else if (!is.null(strip.prefix)) {
    colnames(data) <- str_stripprefix(colnames(data), strip.prefix)
  }
  ## Make sure input and output have the same logginess
  if (input.log && !output.log) {
    data <- 2 ^ data
    input.log <- FALSE
  } else if (!input.log && output.log) {
    data <- log2(data)
    input.log <- TRUE
  }

  if (output.log) {
    fc.operation <- `-`
  } else {
    fc.operation <- `/`
  }

  if (ncol(data) < 2) {
    stop("Need at least 2 columns")
  } else if (ncol(data) == 2) {
    res <- data.frame(row.names=rownames(data), FC=fc.operation(data[,2],data[,1]))
  } else {
    res <- data.frame(llply(colnames(data)[-1], function(col) {
      fc.operation(data[,col], data[,1])
    }))
    names(res) <- str_c("FC.", colnames(data)[-1])
  }
  if (output.log) {
    names(res) <- str_c("log", names(res))
  }
  res
}

do.bayseq.analysis <- function(expdata, feature.counts,
                               replicatecol, groupcols,
                               feature.annot=NULL, splitcols=NULL,
                               sample.subset=NULL, feature.subset=NULL,
                               libsizefactors=getLibsizes,
                               cluster=NULL) {
  force(cluster)
  ## Check dimension matches
  stopifnot(nrow(expdata) == ncol(feature.counts))
  stopifnot(is.null(feature.annot) || nrow(feature.counts) == nrow(feature.annot))
  groupcols <- as.character(groupcols)
  stopifnot(length(groupcols) > 0)

  ## Take appropriate subsets of data and metadata
  if (is.function(sample.subset)) {
    sample.subset <- sample.subset(expdata)
  }
  if (!is.null(sample.subset)) {
    tsmsg("Subsetting samples")
    orig.sample.count <- nrow(expdata)
    expdata <- expdata[sample.subset,]
    feature.counts <- feature.counts[,sample.subset]
    if (is.vector(libsizefactors)) {
      libsizefactors <- libsizefactors[sample.subset]
    }
    tsmsgf("Selected %s out of %s samples.", nrow(expdata), orig.sample.count)
  }

  if (is.function(feature.subset)) {
    feature.subset <- feature.subset(feature.counts, feature.annot)
  }
  if (!is.null(feature.subset)) {
    tsmsg("Subsetting features.")
    orig.feature.count <- nrow(feature.counts)
    feature.counts <- feature.counts[feature.subset,]
    feature.annot <- feature.annot[feature.subset,,drop=FALSE]
    tsmsgf("Selected %s out of %s features.", nrow(feature.counts), orig.feature.count)
  }

  ## If splitcols was provided, split and call recursively on each
  ## sample set
  if (!is.null(splitcols)) {
    tsmsgf("Splitting samples on expdata columns:\n%s", sprint(colnames(expdata[splitcols])))
    splitfactor <- interaction(expdata[splitcols], drop=TRUE)
    splitselections <- split(seq(nrow(expdata)), splitfactor)

    llply(names(splitselections), function(selname) {
      tsmsgf("Doing baySeq analysis on split group: %s", selname)
      sel <- splitselections[[selname]]
      if (is.vector(libsizefactors)) {
        libsizefactors <- libsizefactors[sel]
      }
      do.bayseq.analysis(expdata=expdata[sel,],
                         feature.counts=feature.counts[,sel],
                         feature.annot=feature.annot,
                         replicatecol=replicatecol,
                         groupcols=groupcols,
                         splitcols=NULL,
                         sample.subset=NULL,
                         feature.subset=NULL,
                         libsizefactors=libsizefactors,
                         cluster=cluster)
    })
  } else {
    within(list(), {
      tsmsgf("Beginning baySeq analysis on %s features in %s samples",
             nrow(feature.counts), ncol(feature.counts))

      ## Capture variables in return list
      expdata <- droplevels(expdata)
      feature.counts <- feature.counts
      feature.annot <- feature.annot
      replicatecol <- replicatecol
      groupcols <- groupcols
      names(groupcols) <- groupcols
      tsmsgf("Generating hypotheses from groups:\n%s", sprint(unname(groupcols)))
      groupings <- generate.all.groups(expdata[groupcols])
      group.types <- with(list(gv=str_replace_all(names(groupings), ":.*$", "")),
                          factor(gv, levels=unique(gv)))
      CD <- new("countData",
                data=feature.counts, replicates=expdata[[replicatecol]],
                groups=groupings, annotation=feature.annot)

      if (is.null(libsizefactors)) {
        libsizefactors <- getLibsizes
      }
      if (is.function(libsizefactors)) {
        tsmsg("Estimating library size scaling factors")
        CD@libsizes <- libsizefactors(CD)
      } else {
        tsmsg("Using provided library size factors")
        CD@libsizes <- libsizefactors
      }

      ## It is CPM instead of FPKM because no gene lengths are
      ## provided. This is consistent with DESeq and edgeR
      tsmsg("Calculating CPM and group fold change values")
      cpm <- calcFPKM.CD(CD)
      colnames(cpm) <- str_c("CPM.", colnames(cpm))
      group.cpm <- llply(expdata[groupcols], calcFPKM.byGroup, cD=CD, .parallel=TRUE)
      for (i in 1:length(group.cpm)) {
        colnames(group.cpm[[i]]) <- str_c("CPM.", colnames(group.cpm[[i]]))
      }
      group.fc <- llply(group.cpm, calcFC, strip.prefix="CPM.")

      tsmsg("Getting priors")
      CD <- getPriors.NB(CD, cl = cluster)
      tsmsgf("Finished estimating %s priors", CD@priorType)
      tsmsgf("Getting posterior likelihoods")
      CD <- getLikelihoods.NB(CD, cl = cluster)

      ## Topcounts by group
      tsmsgf("Choosing top DE genes for each group")
      allcounts <- topCounts.by.group.type(CD, group.types, number=nrow(CD), normaliseData=TRUE)[levels(group.types) != "NDE"]
      tsmsg("Reformatting results")
      results <- llply(groupcols, function(group) {
        x <- allcounts[[group]]
        ## Remove per-sample columns
        x <- x[!names(x) %in% colnames(CD@data)]
        ## Add per-group columns and inter-group fold-change columns
        ## (both log and nonlog)
        group.and.fc.data <- data.frame(group.cpm[[group]][rownames(x),],
                                        group.fc[[group]][rownames(x),,drop=FALSE])
        group.and.fc.logdata <- llply(group.and.fc.data, log2)
        names(group.and.fc.logdata) <- str_c("log", names(group.and.fc.data))
        x <- data.frame(group.and.fc.data, x, group.and.fc.logdata)
        ## Put FDR at the front
        x <- reorder.columns(x, start=c("Likelihood", "FDR"))
        ## Put rows back in original order instead of sorting by FDR
        x[row.names(CD@data),]
      })
      topresults <- llply(results, function(x) {
        x <- x[x$FDR <= 0.1,]
        x[order(x$FDR),]
      })

      tsmsgf("Completed baySeq analysis on %s features in %s samples",
             nrow(feature.counts), ncol(feature.counts))
    })
  }
}

formula.from.vector <- function(x, type="deseq") {
  formula.template <- c(deseq="count~%s",
                        edger="~%s",
                        limma="~%s")[type]
  as.formula(sprintf(formula.template, str_c(as.character(x), collapse="+")))
}

as.deseq.formula <- function(f) {
  flist <- as.list(f)
  rhs <- deparse(flist[[length(flist)]])
  newf <- as.formula(sprintf("count ~ %s", rhs))
  environment(newf) <- environment(f)
  newf
}

as.edger.formula <- function(f) {
  flist <- as.list(f)
  rhs <- deparse(flist[[length(flist)]])
  newf <- as.formula(sprintf("~ %s", rhs))
  environment(newf) <- environment(f)
  newf
}

rhs.variables.from.deseq.formula <- function(f) {
  setdiff(all.vars(f), "count")
}

term.labels <- function(f) {
  attr(terms(f), "term.labels")
}

get.group.factor <- function(terms, expdata) {
  if (length(terms) > 0) {
    term.factors <- llply(terms, function(term) as.factor(eval(parse(text=term), envir=expdata)))
    interaction(term.factors, sep=".")
  } else {
    ## With no terms, there is just one big group
    factor(rep(1, nrow(expdata)))
  }
}

## Use lower.tail option instead of subtracting from 1. Maybe gives
## better precision.
nbinomGLMTest <- function (resFull, resReduced) {
  pchisq(resReduced$deviance - resFull$deviance,
         attr(resReduced, "df.residual") - attr(resFull, "df.residual"),
         lower.tail=FALSE,
         log.p=FALSE)
}

quantify.by.group.CountDataSet <- function(cds, group=pData(cds)[[ncol(pData(cds))]], .parallel=TRUE) {
  groupcols <- split(1:ncol(cds), droplevels(as.factor(group)))
  res <- data.frame(llply(groupcols, function(x) {
    subcounts <- counts(cds)[,x]
    subsizefac <- sizeFactors(cds)[x]
    getBaseMeansAndVariances(subcounts, subsizefac)$baseMean
  }, .parallel=.parallel),
             row.names=row.names(counts(cds)))
  names(res) <- str_c("CPM.", names(res))
  res
}

plotDispersionTrend.CountDataSet <- function(cds) {
  x <- data.frame(CPM=getBaseMeansAndVariances(counts(cds), sizeFactors(cds))$baseMean,
                  fitInfo(cds)[c("fittedDispEsts", "perGeneDispEsts")])
  ggplot(x) + aes(x=CPM) + geom_hex(aes(y=perGeneDispEsts)) +
    geom_line(aes(y=fittedDispEsts)) +
      scale_x_log10() + scale_y_log10() + ylab("Dispersion") + xlab("CPM") + ggtitle("Dispersion Trend")
}

plotDispersionTrend.DGEList <- function(dge) {
  x <- data.frame(#CPM=getBaseMeansAndVariances(dge$counts, dge$samples$norm.factors)$baseMean,
                  CPM=2 ^ dge$logCPM,
                  trend=dge$trended.dispersion,
                  tag=dge$tagwise.dispersion)
  ggplot(x) + aes(x=CPM) + geom_hex(aes(y=tag)) +
    geom_line(aes(y=trend)) +
      scale_x_log10() + scale_y_log10() + ylab("Dispersion") + xlab("CPM") + ggtitle("Dispersion Trend")
}

## Allow df.residual to be stored in attribute (to support DESeq)
df.residual.default <- function(object, ...) {
  if ("df.residual" %in% names(attributes(object))) {
    attr(object, "df.residual")
  } else {
    stats:::df.residual.default(object, ...)
  }
}

## Make gof generic for both edgeR and DESeq
gof <- function(object, pcutoff = 0.1, adjust = "holm", plot = FALSE,
                main = "qq-plot of genewise goodness of fit", ...) {
  UseMethod("gof")
}
gof.default <- function(object, pcutoff = 0.1, adjust = "holm", plot = FALSE,
                main = "qq-plot of genewise goodness of fit", ...) {
  if (is.numeric(object)) {
    gof.stats <- object
    gof.names <- names(object)
  } else {
    gof.stats <- deviance(object)
    gof.names <- row.names(object)
  }
  gof.pvals <- pchisq(gof.stats, df = df.residual(object), lower.tail = FALSE,
                      log.p = FALSE)
  gof.padj <- p.adjust(gof.pvals, method = adjust)
  outlier <- gof.padj < pcutoff
  if (plot) {
    n <- length(gof.stats)
    z <- zscoreGamma(gof.stats, shape = df.residual(object)/2,
                     scale = 2)
    col <- rep("black", n)
    col[outlier] <- "blue"
    pch <- rep(1, n)
    pch[outlier] <- 16
    data <- data.frame(x=qnorm(ppoints(n))[order(order(z))],
                       y=z,
                       p=gof.pvals,
                       q=p.adjust(gof.pvals, method = adjust),
                       outlier=outlier,
                       col=col,
                       pch=pch)
    print(ggplot(data) + aes(x=x, y=y, color=outlier, shape=outlier) +
          geom_point() +
          scale_colour_manual(values=c("blue", "red")) +
          scale_shape() +
          geom_abline(slope=1, intercept=0))
    ## print(ggplot(data) +
    ##       aes(x=x, y=y, color=outlier, shape=outlier) +
    ##       geom_point() + geom_abline(slope=1, intercept=0) +
    ##       scale_colour_manual(values=c("blue", "red")))
    ## qqnorm(z, col = col, pch = pch, main = main, ...)
    ## abline(0, 1)
  }
  data.frame(row.names=make.unique(gof.names), deviance = gof.stats,
             pvalue = gof.pvals, padj = gof.padj,
             outlier = outlier, df = df.residual(object)[1])
}

do.deseq.analysis <- function(expdata, feature.counts,
                              null.formula=~1, alt.formula=~condition,
                              disp.formula=alt.formula,
                              feature.annot=NULL, splitcols=NULL,
                              sample.subset=NULL, feature.subset=NULL,
                              libsizefactors=estimateSizeFactors,
                              .parallel=TRUE, ...) {
  ## Check dimension matches
  stopifnot(nrow(expdata) == ncol(feature.counts))
  stopifnot(is.null(feature.annot) || nrow(feature.counts) == nrow(feature.annot))

  ## Convert vectors of variables to real formulas. i.e. convert
  ## c("a", "b", "c") into count~a+b+c.
  if (is.vector(null.formula)) {
    null.formula <- formula.from.vector(null.formula, "deseq")
  }
  if (is.vector(alt.formula)) {
    alt.formula <- formula.from.vector(alt.formula, "deseq")
  }
  if (is.vector(disp.formula)) {
    disp.formula <- formula.from.vector(disp.formula, "deseq")
  } else if (is.null(disp.formula)) {
    disp.formula <- alt.formula
  }

  ## Take appropriate subsets of data and metadata
  if (is.function(sample.subset)) {
    sample.subset <- sample.subset(expdata)
  }
  if (!is.null(sample.subset)) {
    tsmsg("Subsetting samples")
    orig.sample.count <- nrow(expdata)
    expdata <- expdata[sample.subset,]
    feature.counts <- feature.counts[,sample.subset]
    if (is.vector(libsizefactors)) {
      libsizefactors <- libsizefactors[sample.subset]
    }
    tsmsgf("Selected %s out of %s samples.", nrow(expdata), orig.sample.count)
  }

  if (is.function(feature.subset)) {
    feature.subset <- feature.subset(feature.counts, feature.annot)
  }
  if (!is.null(feature.subset)) {
    tsmsg("Subsetting features.")
    orig.feature.count <- nrow(feature.counts)
    feature.counts <- feature.counts[feature.subset,]
    feature.annot <- feature.annot[feature.subset,,drop=FALSE]
    tsmsgf("Selected %s out of %s features.", nrow(feature.counts), orig.feature.count)
  }

  ## If splitcols was provided, split and call recursively on each
  ## sample set
  if (!is.null(splitcols)) {
    tsmsgf("Splitting samples on expdata columns:\n%s", sprint(colnames(expdata[splitcols])))
    splitfactor <- interaction(expdata[splitcols], drop=TRUE)
    splitselections <- split(seq(nrow(expdata)), splitfactor)

    llply(names(splitselections), function(selname) {
      tsmsgf("Doing DESeq analysis on split group: %s", selname)
      sel <- splitselections[[selname]]
      if (is.vector(libsizefactors)) {
        libsizefactors <- libsizefactors[sel]
      }
      do.deseq.analysis(expdata=expdata[sel,],
                        feature.counts=feature.counts[,sel],
                        feature.annot=feature.annot,
                        null.formula=null.formula,
                        alt.formula=alt.formula,
                        disp.formula=disp.formula,
                        splitcols=NULL,
                        sample.subset=NULL,
                        feature.subset=NULL,
                        libsizefactors=libsizefactors,
                        .parallel=FALSE, ...)
    }, .parallel=.parallel)
  } else {
    extra.args <- list(...)

    within(list(), {
      tsmsgf("Beginning DESeq analysis on %s features in %s samples",
             nrow(feature.counts), ncol(feature.counts))

      ## Capture variables in return list
      expdata <- droplevels(expdata)
      feature.counts <- feature.counts
      feature.annot <- feature.annot
      null.formula <- null.formula
      alt.formula <- alt.formula
      disp.formula <- disp.formula
      tsmsgf("Testing hypotheses:\n\tNull: %s\n\tAlternate: %s",
             deparse(null.formula), deparse(alt.formula))
      used.variables <- union(rhs.variables.from.deseq.formula(null.formula),
                              rhs.variables.from.deseq.formula(alt.formula))
      alt.terms <- setdiff(term.labels(alt.formula),
                           term.labels(null.formula))
      disp.dropped.terms <- setdiff(term.labels(alt.formula),
                                    term.labels(disp.formula))
      groupvar <- get.group.factor(alt.terms, expdata)
      stopifnot(nlevels(groupvar) > 1)

      design.null <- model.matrix(null.formula, data=expdata)
      design.alt <- model.matrix(alt.formula, data=expdata)
      design.disp <- model.matrix(disp.formula, data=expdata)
      stopifnot(all(colnames(design.null) %in% colnames(design.alt)))
      stopifnot(all(colnames(design.disp) %in% colnames(design.alt)))
      design.alt.cols <- setdiff(colnames(design.alt),
                                 colnames(design.null))

      cds <- newCountDataSet(feature.counts, expdata[used.variables],
                             featureData=if (!is.null(feature.annot)) as(feature.annot, "AnnotatedDataFrame"))
      ## Map(factor, pData(cds)[used.variables])

      if (is.null(libsizefactors)) {
        libsizefactors <- estimateSizeFactors
      }
      if (is.function(libsizefactors)) {
        tsmsg("Estimating library size scaling factors")
        libsize.results <- libsizefactors(cds)
        ## Handle returning either the whole cds object or just a vector of libsizefactors
        if (is(libsize.results, class(cds))) {
          cds <- libsize.results
        } else if (is.vector(libsize.results)) {
          sizeFactors(cds) <- libsize.results
        } else {
          stop("Library size function returned something unknown")
        }
      } else {
        tsmsg("Using provided library size factors")
        sizeFactors(cds) <- libsizefactors
      }

      default.dispersion.args <- list(method="pooled-CR", fitType="parametric",
                                      modelFormula=disp.formula, modelFrame=expdata[used.variables])
      dispersion.args <- extra.args[str_detect(names(extra.args), "^disp\\.")]
      names(dispersion.args) <- str_replace(names(dispersion.args), "^disp\\.", "")
      dispersion.args <- merge.defaults.into.list(dispersion.args,
                                                  defaults=default.dispersion.args)
      dispersion.args$modelFormula <- as.deseq.formula(dispersion.args$modelFormula)

      tsmsgf("Running \"estimateDispersions\" with arguments:\n%s", sprint(dispersion.args))
      if (length(disp.dropped.terms) > 0) {
        tsmsgf("Dropping terms %s for dispersion estimation", deparse(unname(disp.dropped.terms)))
      }
      cds <- do.call(estimateDispersions,
                     c(list(object=cds), dispersion.args))

      tsmsg("Fitting GLMs")
      fit <- glmFit.CountDataSet(cds, design.alt)
      gof.info <- gof(fit)

      tsmsg("Performing likelihood ratio tests")
      lrt <- glmLRT(fit, coef=design.alt.cols)

      tsmsg("Creating result tables")
      results <- data.frame(topTags(lrt, n=nrow(lrt)))[featureNames(cds),]

      ## Eval this in a temp list so we don't capture these
      ## variables
      temp <- within(list(), {
        tsmsg("Calculating means and fold changes")

        ## Calculate overall mean and the mean for each group
        meanCPM <- 2 ^ results$logCPM
        groupCPM <- quantify.by.group.CountDataSet(cds, group=groupvar)
        grouplogCPM <- data.frame(log2(as.matrix(groupCPM)))
        names(grouplogCPM) <- str_c("log", names(groupCPM))

        ## Get fold changes and log fold changes
        logFC <- results[str_detect(names(results), "^logFC")]
        FC <- data.frame(2 ^ as.matrix(logFC))
        names(FC) <- str_replace(names(logFC), "^logFC", "FC")

        res <- data.frame(row.names=row.names(feature.annot),
                          CPM.mean=meanCPM,
                          groupCPM, FC)
        logres <- data.frame(row.names=row.names(res),
                             logCPM.mean=results$logCPM,
                             grouplogCPM, logFC)
      })

      results <- data.frame(temp$res,
                            GOF.outlier=gof.info$outlier,
                            results[setdiff(names(results), c("logCPM", names(temp$res), names(temp$logres)))],
                            temp$logres)
      results <- reorder.columns(results, start=c("LR", "PValue", "FDR"))

      tsmsg("Choosing top DE genes")
      topresults <- results[!is.na(results$FDR) & results$FDR <= 0.1,]
      topresults <- topresults[order(topresults$FDR,
                                     topresults$PValue,
                                     -topresults$LR),]

      tsmsgf("Completed DESeq analysis on %s features in %s samples",
             nrow(feature.counts), ncol(feature.counts))
    })
  }
}

## EdgeR
makeEdgeRDispersionEstimator <- function(trended=TRUE, tagwise=TRUE, glm=TRUE,
                                         common.opts=list(),
                                         trend.opts=list(),
                                         tagwise.opts=list(),
                                         return.type=c("DGEList", "vector")) {
  return.type <- match.arg(return.type)
  if (glm) {
    common.dispfun <- curry(estimateGLMCommonDisp, common.opts)
    trended.dispfun <- curry(estimateGLMTrendedDisp, trend.opts)
    tagwise.dispfun <- curry(estimateGLMTagwiseDisp, tagwise.opts)
    function(dge, design, offset=NULL) {
      if (trended) {
        dge <- trended.dispfun(dge, design, offset)
      } else {
        dge <- common.dispfun(dge, design, offset)
      }
      if (tagwise) {
        dge <- tagwise.dispfun(dge, design, offset)
      }
      if (return.type == "DGEList") {
        dge
      } else if (return.type == "vector") {
        if (tagwise) {
          dge$tagwise.dispersion
        } else if (trended) {
          dge$trended.dispersion
        } else {
          dge$common.dispersion
        }
      }
    }
  } else {
    common.dispfun <- curry(estimateCommonDisp, common.opts)
    trended.dispfun <- curry(estimateTrendedDisp, trend.opts)
    tagwise.dispfun <- curry(estimateTagwiseDisp, tagwise.opts)
    function(dge) {
      if (trended) {
        dge <- trended.dispfun(dge)
      } else {
        dge <- common.dispfun(dge)
      }
      if (tagwise) {
        dge <- tagwise.dispfun(dge)
      }
      if (return.type == "DGEList") {
        dge
      } else if (return.type == "vector") {
        if (tagwise) {
          dge$tagwise.dispersion
        } else if (trended) {
          dge$trended.dispersion
        } else {
          dge$common.dispersion
        }
      }
    }
  }
}

setMethod("estimateDispersions", signature=(object="DGEList"),
          function(object, ...) makeEdgeRDispersionEstimator()(dge=object, ...))

run.aqm <- function(expdata, feature.counts,
                    feature.subset=NULL, sample.subset=NULL,
                    reporttitle="arrayQualityMetrics report for ",
                    outdir=reporttitle,
                    splitcols=NULL,
                    reporttitle.template="arrayQualityMetrics report for %s",
                    outdir.template=reporttitle.template,
                    .parallel=TRUE, ...) {
  ## Take appropriate subsets of data and metadata
  if (is.function(sample.subset)) {
    sample.subset <- sample.subset(expdata)
  }
  if (!is.null(sample.subset)) {
    expdata <- expdata[sample.subset,]
    feature.counts <- feature.counts[,sample.subset]
  }

  if (is.function(feature.subset)) {
    feature.subset <- feature.subset(feature.counts, feature.annot)
  }
  if (!is.null(feature.subset)) {
    feature.counts <- feature.counts[feature.subset,]
    feature.annot <- feature.annot[feature.subset,,drop=FALSE]
  }

  ## If splitcols was provided, split and call recursively on each
  ## sample set
  if (!is.null(splitcols)) {
    splitfactor <- interaction(expdata[splitcols], drop=TRUE)
    splitselections <- split(seq(nrow(expdata)), splitfactor)

    llply(names(splitselections), function(selname) {
      sel <- splitselections[selname]
      sel.title <- sprintf(reporttitle.template, selname)
      sel.outdir <- sprintf(outdir.template, selname)
      run.aqm(expdata=expdata[sel,],
              feature.counts=feature.counts[,sel],
              reporttitle=seltitle, outdir=sel.outdir,
              .parallel=FALSE, ...)
    }, .parallel=.parallel)
  } else {
    cds <- newCountDataSet(feature.counts, expdata)
    eset <- as(cds, "ExpressionSet")
    if (!file.exists(outdir)) {
      dir.create(outdir, recursive=TRUE)
    }
    arrayQualityMetrics(eset, reporttitle=reporttitle, outdir=outdir, ...)
  }
}

## nbinomGLMTest work-alike for edgeR DGEGLM class
nbinomGLMTest.DGEGLM <- function (resFull, resReduced) {
  pchisq(resReduced$deviance - resFull$deviance,
         resReduced$df.residual - resFull$df.residual,
         lower.tail=FALSE,
         log.p=FALSE)
}

quantify.by.group.DGEList <- function(dge, group=dge$samples$group, design=NULL,
                                      offset=NULL,
                                      dispersion.estimator=makeEdgeRDispersionEstimator(),
                                      .parallel=TRUE) {
  groupcols <- split(1:ncol(dge), droplevels(group))
  names(groupcols) <- str_c("CPM.", names(groupcols))
  data.frame(llply(groupcols, function(cols) {
    subdge <- DGEList(counts=dge$counts[,cols],
                      genes=dge$genes,
                      group=dge$samples$group[cols],
                      lib.size=dge$samples$lib.size[cols],
                      norm.factors=dge$samples$norm.factors[cols])
    if (is.null(design)) {
      subdesign <- NULL
    } else {
      subdesign <- design[cols,]
    }
    subdesign <- cbind(rep(1, length(cols)))
    subdge <- dispersion.estimator(subdge, subdesign, offset)
    subfit <- glmFit(subdge, subdesign)
    res <- (subfit$abundance + log(1e+06)) / log(2)
    ## EdgeR gives log2 values by default, so exponentiate
    2^res
  }, .parallel=.parallel),
             row.names=row.names(dge))
}

## This version of glmQLFTest excludes genes with zero counts in all
## samples (logCPM == -Inf), because these cause squeezeVar to throw
## errors. Instead, it just arbitrarily sets their F statistics to
## zero and P-Values to 1, which is appropriate because obviously we
## cannot call differential expression if all counts are zero.
glmQLFTest.safe <- function (glmfit, coef = ncol(glmfit$design), contrast = NULL,
                             abundance.trend = TRUE)
{
  present <- is.finite(glmfit$abundance) & is.finite(glmfit$dispersion)
  out <- glmLRT(glmfit, coef = coef, contrast = contrast)
  df.residual <- glmfit$df.residual
  s2 <- glmfit$deviance/df.residual
  df.residual[s2 < 1e-14] <- 0
  s2 <- pmax(s2, 0)
  if (abundance.trend)  {
    s2.fit <- squeezeVar(s2[present],
                         df = df.residual[present],
                         covariate = glmfit$abundance[present])
  }
  else {
    s2.fit <- squeezeVar(s2, df = df.residual)
  }
  F <- rep(0, length(s2))
  F[present] <- (out$table$LR/out$df.test)[present]/s2.fit$var.post
  df.total <- s2.fit$df.prior + df.residual
  max.df.residual <- ncol(glmfit$counts) - ncol(glmfit$design)
  df.total <- min(df.total, length(s2) * max.df.residual)
  F.pvalue <- pf(F, df1 = out$df.test, df2 = df.total, lower.tail = FALSE,
                 log.p = FALSE)
  i <- s2.fit$var.post < 1
  if (any(i)) {
    chisq.pvalue <- pchisq(out$table$LR[i], df = out$df.test[i],
                           lower.tail = FALSE, log.p = FALSE)
    F.pvalue[i] <- pmax(F.pvalue[i], chisq.pvalue)
  }
  out$table$LR <- out$table$PValue <- NULL
  out$table$F <- F
  out$table$PValue <- F.pvalue
  out$s2.fit <- s2.fit
  out$df.prior <- s2.fit$df.prior
  out$df.total <- df.total
  out
}

do.edger.analysis <- function(expdata, feature.counts,
                              null.formula=~1, alt.formula=~condition,
                              disp.formula=alt.formula,
                              feature.annot=NULL, splitcols=NULL,
                              sample.subset=NULL, feature.subset=NULL,
                              dispersion.estimator=makeEdgeRDispersionEstimator(),
                              libsizefactors="TMM",
                              test.method=c("QL", "LR"),
                              offset=NULL,
                              .parallel=TRUE) {
  test.method <- match.arg(test.method)

  ## Check dimension matches
  stopifnot(nrow(expdata) == ncol(feature.counts))
  stopifnot(is.null(feature.annot) || nrow(feature.counts) == nrow(feature.annot))

  ## Allow specifying any valid arg to calcNormFactors(method=...)
  if (is.character(libsizefactors)) {
    if (libsizefactors %in% eval(formals(calcNormFactors)$method)) {
      tsmsgf("Using method %s to calculate library normalization factors", libsizefactors)
      libsizefactors <- curry(calcNormFactors, list(method=libsizefactors[[1]]))
    } else {
      stop(sprintf("Unknown library size estimation method: %s", libsizefactors[[1]]))
    }
  }

  ## Expand offset to same dimensions as counts
  if (!is.null(offset)) {
    tsmsgf("Using provided offsets in GLM fits and dispersion estimation.")
    offset <- expandAsMatrix(offset, dim(feature.counts))
  }

  ## Convert vectors of variables to real formulas. i.e. convert
  ## c("a", "b", "c") into count~a+b+c.
  if (is.vector(null.formula)) {
    null.formula <- formula.from.vector(null.formula, "edger")
  }
  if (is.vector(alt.formula)) {
    alt.formula <- formula.from.vector(alt.formula, "edger")
  }
  if (is.vector(disp.formula)) {
    disp.formula <- formula.from.vector(disp.formula, "deseq")
  } else if (is.null(disp.formula)) {
    disp.formula <- alt.formula
  }

  ## Take appropriate subsets of data and metadata
  if (is.function(sample.subset)) {
    sample.subset <- sample.subset(expdata)
  }
  if (!is.null(sample.subset)) {
    tsmsg("Subsetting samples")
    orig.sample.count <- nrow(expdata)
    expdata <- expdata[sample.subset,]
    feature.counts <- feature.counts[,sample.subset]
    offset <- offset[,sample.subset]
    if (is.vector(libsizefactors)) {
      libsizefactors <- libsizefactors[sample.subset]
    }
    tsmsgf("Selected %s out of %s samples.", nrow(expdata), orig.sample.count)
  }

  if (is.function(feature.subset)) {
    feature.subset <- feature.subset(feature.counts, feature.annot)
  }
  if (!is.null(feature.subset)) {
    tsmsg("Subsetting features.")
    orig.feature.count <- nrow(feature.counts)
    feature.counts <- feature.counts[feature.subset,]
    offset <- offset[feature.subset,]
    feature.annot <- feature.annot[feature.subset,,drop=FALSE]
    tsmsgf("Selected %s out of %s features.", nrow(feature.counts), orig.feature.count)
  }

  ## If splitcols was provided, split and call recursively on each
  ## sample set
  if (!is.null(splitcols)) {
    tsmsgf("Splitting samples on expdata columns:\n%s", sprint(colnames(expdata[splitcols])))
    splitfactor <- interaction(expdata[splitcols], drop=TRUE)
    splitselections <- split(seq(nrow(expdata)), splitfactor)

    llply(names(splitselections), function(selname) {
      tsmsgf("Doing edgeR analysis on split group: %s", selname)
      sel <- splitselections[[selname]]
      if (is.vector(libsizefactors)) {
        libsizefactors <- libsizefactors[sel]
      }
      do.edger.analysis(expdata=expdata[sel,],
                        feature.counts=feature.counts[,sel],
                        offset=offset[,sel],
                        feature.annot=feature.annot,
                        null.formula=null.formula,
                        alt.formula=alt.formula,
                        disp.formula=disp.formula,
                        splitcols=NULL,
                        sample.subset=NULL,
                        feature.subset=NULL,
                        libsizefactors=libsizefactors,
                        dispersion.estimator=dispersion.estimator,
                        test.method=test.method,
                        .parallel=FALSE)
    }, .parallel=.parallel)
  } else {
    within(list(), {
      tsmsgf("Beginning edgeR analysis on %s features in %s samples",
             nrow(feature.counts), ncol(feature.counts))

      ## Capture variables in return list
      expdata <- droplevels(expdata)
      feature.counts <- feature.counts
      offset <- offset
      feature.annot <- feature.annot
      null.formula <- null.formula
      alt.formula <- alt.formula
      disp.formula <- disp.formula
      tsmsgf("Testing hypotheses:\n\tNull: %s\n\tAlternate: %s",
             deparse(null.formula), deparse(alt.formula))
      used.variables <- union(rhs.variables.from.deseq.formula(null.formula),
                              rhs.variables.from.deseq.formula(alt.formula))
      alt.terms <- setdiff(term.labels(alt.formula),
                           term.labels(null.formula))
      disp.dropped.terms <- setdiff(term.labels(alt.formula),
                                    term.labels(disp.formula))
      groupvar <- get.group.factor(alt.terms, expdata)

      design.null <- model.matrix(null.formula, data=expdata)
      design.alt <- model.matrix(alt.formula, data=expdata)
      design.disp <- model.matrix(disp.formula, data=expdata)
      stopifnot(all(colnames(design.null) %in% colnames(design.alt)))
      stopifnot(all(colnames(design.disp) %in% colnames(design.alt)))
      design.alt.cols <- setdiff(colnames(design.alt),
                                 colnames(design.null))

      dge <- DGEList(counts=feature.counts[,rownames(expdata)],
                     group=groupvar,
                     genes=feature.annot)
      if (is.null(libsizefactors)) {
        libsizefactors <- calcNormFactors
      }
      if (is.function(libsizefactors)) {
        tsmsg("Estimating library size scaling factors")
        libsize.results <- libsizefactors(dge)
        ## Handle returning either a modified object or just a vector
        ## of size factors
        if (is(libsize.results, class(dge))) {
          dge <- libsize.results
        } else if (is.vector(libsize.results)) {
          dge$samples$norm.factors <- libsize.results
        } else {
          stop("Library size function returned something unknown")
        }
      } else {
        tsmsg("Using provided library size factors")
        dge$samples$norm.factors <- libsizefactors
      }
      tsmsg("Estimating dispersions")
      if (length(disp.dropped.terms) > 0) {
        tsmsgf("Dropping terms %s for dispersion estimation", deparse(unname(disp.dropped.terms)))
      }
      dge <- dispersion.estimator(dge, design.disp, offset)

      tsmsg("Fitting GLMs")
      fit <- glmFit(dge, design.alt, offset=offset)
      gof.info <- gof(fit)

      tsmsg("Performing likelihood ratio tests")
      lrt <- glmLRT(fit, coef=design.alt.cols)
      tsmsg("Performing quasi-likelihood F-tests")
      qlft <- glmQLFTest.safe(fit, coef=design.alt.cols)

      tsmsg("Creating result tables")
      results.lr <- data.frame(topTags(lrt, n=nrow(lrt)))[rownames(dge),]
      results.ql <- data.frame(topTags(qlft, n=nrow(qlft)))[rownames(dge),]

      ## Eval this in a temp list so we don't capture these
      ## variables
      temp <- within(list(), {
        tsmsg("Calculating means and fold changes")

        ## Calculate overall mean and the mean for each group
        meanCPM <- 2 ^ results.lr$logCPM
        groupCPM <- quantify.by.group.DGEList(dge, group=groupvar, design=design.alt,
                                              offset=offset,
                                              dispersion.estimator=dispersion.estimator,
                                              .parallel=.parallel)
        grouplogCPM <- data.frame(log2(as.matrix(groupCPM)))
        names(grouplogCPM) <- str_c("log", names(groupCPM))

        ## Get fold changes and log fold changes
        logFC <- results.lr[str_detect(names(results.lr), "^logFC")]
        FC <- data.frame(2 ^ as.matrix(logFC))
        names(FC) <- str_replace(names(logFC), "^logFC", "FC")

        res <- data.frame(row.names=row.names(feature.annot),
                          CPM.mean=meanCPM,
                          groupCPM, FC)
        logres <- data.frame(row.names=row.names(res),
                             logCPM.mean=results.lr$logCPM,
                             grouplogCPM, logFC)
      })

      results.lr <- data.frame(temp$res,
                            GOF.outlier=gof.info$outlier,
                            results.lr[setdiff(names(results.lr), c("logCPM", names(temp$res), names(temp$logres)))],
                            temp$logres)
      results.lr <- reorder.columns(results.lr, start=c("LR", "PValue", "FDR"))

      results.ql <- data.frame(temp$res,
                            GOF.outlier=gof.info$outlier,
                            results.ql[setdiff(names(results.ql), c("logCPM", names(temp$res), names(temp$logres)))],
                            temp$logres)
      results.ql <- reorder.columns(results.ql, start=c("F", "PValue", "FDR"))

      tsmsg("Choosing top DE genes")
      topresults.lr <- results.lr[!is.na(results.lr$FDR) & results.lr$FDR <= 0.1,]
      topresults.lr <- topresults.lr[order(topresults.lr$FDR,
                                           topresults.lr$PValue,
                                           -topresults.lr$LR),]
      topresults.ql <- results.ql[!is.na(results.ql$FDR) & results.ql$FDR <= 0.1,]
      topresults.ql <- topresults.ql[order(topresults.ql$FDR,
                                        topresults.ql$PValue,
                                        -topresults.ql$F),]
      switch(test.method,
             LR={
               tsmsg("Using likelihood ratio test results")
               test <- lrt
               results <- results.lr
               topresults <- topresults.lr
               test.stat <- "LR"
             },
             QL={
               tsmsg("Using quasi-likelihood F-test results")
               test <- qlft
               results <- results.ql
               topresults <- topresults.ql
               test.stat <- "F"
             })

      tsmsgf("Completed edgeR analysis on %s features in %s samples",
             nrow(feature.counts), ncol(feature.counts))
    })
  }
}

## Replacement for eBayes that does unmoderated t-tests. Not
## recommended for actual use. See:
## https://stat.ethz.ch/pipermail/bioconductor/2010-September/035210.html
no.eBayes <- function(fit, ...) {
  warning("Unmoderated t-statistics are not recommended for use. Please use eBayes instead.")

  ## We call eBayes here to populate
  fit <- eBayes(fit, ...)
  fit$df.prior <- 0
  ## fit$s2.prior <- NA
  ## fit$var.prior <- NA
  ## fit$proportion <- NA
  ## fit$s2.post <- NA
  fit$t <- fit$coef / fit$stdev.unscaled / fit$sigma
  fit$df.total <- fit$df.residual
  fit$p.value <- pt(-abs(fit$t), df = fit$df.residual)
  ## fit$lods <- NA
  if (!is.null(fit$design) && is.fullrank(fit$design)) {
    F.stat <- classifyTestsF(fit, fstat.only = TRUE)
    fit$F <- as.vector(F.stat)
    df1 <- attr(F.stat, "df1")
    df2 <- attr(F.stat, "df2")
    if (df2[1] > 1e+06)
      fit$F.p.value <- pchisq(df1 * fit$F, df1, lower.tail = FALSE)
    else fit$F.p.value <- pf(fit$F, df1, df2, lower.tail = FALSE)
  }
  fit
}

do.limma.analysis <- function(expdata, feature.counts,
                              null.formula=~1, alt.formula=~condition, block.cols=NULL,
                              feature.annot=NULL, splitcols=NULL,
                              sample.subset=NULL, feature.subset=NULL,
                              libsizefactors=calcNormFactors,
                              transformation=c("voom", "vst", "none"),
                              do.ebayes=TRUE,
                              .parallel=TRUE) {
  transformation <- match.arg(transformation)

  ## Check dimension matches
  stopifnot(nrow(expdata) == ncol(feature.counts))
  stopifnot(is.null(feature.annot) || nrow(feature.counts) == nrow(feature.annot))

  ## Convert vectors of variables to real formulas. i.e. convert
  ## c("a", "b", "c") into count~a+b+c.
  if (is.vector(null.formula)) {
    null.formula <- formula.from.vector(null.formula, "limma")
  }
  if (is.vector(alt.formula)) {
    alt.formula <- formula.from.vector(alt.formula, "limma")
  }

  ## Take appropriate subsets of data and metadata
  if (is.function(sample.subset)) {
    sample.subset <- sample.subset(expdata)
  }
  if (!is.null(sample.subset)) {
    tsmsg("Subsetting samples")
    orig.sample.count <- nrow(expdata)
    expdata <- expdata[sample.subset,]
    feature.counts <- feature.counts[,sample.subset]
    if (is.vector(libsizefactors)) {
      libsizefactors <- libsizefactors[sample.subset]
    }
    tsmsgf("Selected %s out of %s samples.", nrow(expdata), orig.sample.count)
  }

  if (is.function(feature.subset)) {
    feature.subset <- feature.subset(feature.counts, feature.annot)
  }
  if (!is.null(feature.subset)) {
    tsmsg("Subsetting features.")
    orig.feature.count <- nrow(feature.counts)
    feature.counts <- feature.counts[feature.subset,]
    feature.annot <- feature.annot[feature.subset,,drop=FALSE]
    tsmsgf("Selected %s out of %s features.", nrow(feature.counts), orig.feature.count)
  }

  ## If splitcols was provided, split and call recursively on each
  ## sample set
  if (!is.null(splitcols)) {
    tsmsgf("Splitting samples on expdata columns:\n%s", sprint(colnames(expdata[splitcols])))
    splitfactor <- interaction(expdata[splitcols], drop=TRUE)
    splitselections <- split(seq(nrow(expdata)), splitfactor)

    llply(names(splitselections), function(selname) {
      tsmsgf("Doing Limma analysis on split group: %s", selname)
      sel <- splitselections[[selname]]
      if (is.vector(libsizefactors)) {
        libsizefactors <- libsizefactors[sel]
      }
      do.limma.analysis(expdata=expdata[sel,],
                        feature.counts=feature.counts[,sel],
                        feature.annot=feature.annot,
                        null.formula=null.formula,
                        alt.formula=alt.formula,
                        splitcols=NULL,
                        sample.subset=NULL,
                        feature.subset=NULL,
                        libsizefactors=libsizefactors,
                        transformation=transformation,
                        do.ebayes=do.ebayes,
                        .parallel=FALSE)
    }, .parallel=.parallel)
  } else {
    within(list(), {
      tsmsgf("Beginning Limma analysis on %s features in %s samples",
             nrow(feature.counts), ncol(feature.counts))

      ## Capture variables in return list
      expdata <- droplevels(expdata)
      block.cols <- block.cols
      feature.counts <- feature.counts
      feature.annot <- feature.annot
      null.formula <- null.formula
      alt.formula <- alt.formula
      tsmsgf("Testing hypotheses:\n\tNull: %s\n\tAlternate: %s",
             deparse(null.formula), deparse(alt.formula))
      used.variables <- union(rhs.variables.from.deseq.formula(null.formula),
                              rhs.variables.from.deseq.formula(alt.formula))
      alt.terms <- setdiff(term.labels(alt.formula),
                           term.labels(null.formula))
      groupvar <- get.group.factor(alt.terms, expdata)

      design.null <- model.matrix(null.formula, data=expdata)
      design.alt <- model.matrix(alt.formula, data=expdata)
      stopifnot(all(colnames(design.null) %in% colnames(design.alt)))
      design.alt.cols <- setdiff(colnames(design.alt),
                                 colnames(design.null))

      if (transformation != "none") {
        dge <- DGEList(counts=feature.counts[,rownames(expdata)],
                       group=groupvar,
                       genes=feature.annot)
        if (is.null(libsizefactors)) {
          libsizefactors <- calcNormFactors
        }
        if (is.function(libsizefactors)) {
          tsmsg("Estimating library size scaling factors")
          libsize.results <- libsizefactors(dge)
          ## Handle returning either a modified object or just a vector
          ## of size factors
          if (is(libsize.results, class(dge))) {
            dge <- libsize.results
          } else if (is.vector(libsize.results)) {
            dge$samples$norm.factors <- libsize.results
          } else {
            stop("Library size function returned something unknown")
          }
        } else {
          tsmsg("Using provided library size factors")
          dge$samples$norm.factors <- libsizefactors
        }
      }
      eset <-
        switch(transformation,
               voom={
                 tsmsg("Transforming counts to weighted log2(CPM) values using limma voom")
                 voom(dge, design.alt, plot=FALSE)
               },
               vst={
                 tsmsg("Using DESeq variance-stabilizing transformation")
                 as(newCountDataSet(feature.counts, expdata[used.variables],
                                    featureData=as(feature.annot, "AnnotatedDataFrame")),
                    "ExpressionSet")
               },
               none={
                 tsmsg("Skipping transformation step")
                 ExpressionSet(feature.counts,
                               phenoData=as(expdata, "AnnotatedDataFrame"),
                               featureData=as(feature.annot, "AnnotatedDataFrame"))
               })

      if (is.null(block.cols)) {
        block.factor <- NULL
        corfit <- NULL
        tsmsg("Fitting linear models")
        fit <- lmFit(eset, design.alt)
      } else {
        tsmsgf("Estimating correlation for blocking factors %s", deparse(unname(block.cols)))
        block.factor <- interaction(expdata[block.cols], drop=TRUE)
        corfit <- duplicateCorrelation(eset, design.alt, block=block.factor)
        tsmsg("Fitting linear models")
        fit <- lmFit(eset, design.alt, block=block.factor, correlation=corfit$consensus.correlation)
      }
      if (do.ebayes) {
        tsmsg("Computing moderated test statistics by eBayes shrinkage")
        fit <- eBayes(fit)
      } else {
        tsmsg("Computing unmoderated test statistics.")
        fit <- no.eBayes(fit)
      }
      tsmsg("Creating result table")
      results <- data.frame(topTable(fit, coef=design.alt.cols,
                                     n=nrow(eset), sort.by="none"),
                            row.names=row.names(fit$genes))
      ## Change names to be consistent with other analysis pipelines
      results <- rename(results, c(adj.P.Val="FDR", P.Value="PValue", AveExpr="logAveExpr"))

      ## Eval this in a temp list so we don't capture these
      ## variables
      temp <- within(list(), {
        tsmsg("Calculating means and fold changes")
        logFC <- results[names(results) %in% c("logFC", colnames(fit$coefficients))]
        ## Ensure each name begins with "logFC"
        names(logFC) <- str_replace(names(logFC),
                                    perl("^(?!logFC(.|$))"),
                                    "logFC.")
        FC <- data.frame(2 ^ as.matrix(logFC))
        names(FC) <- str_replace(names(logFC), "^logFC", "FC")

        res <- data.frame(row.names=row.names(feature.annot),
                          AveExpr=2^results$AveExpr,
                          FC)
        logres <- data.frame(row.names=row.names(res),
                             logAveExpr=results$AveExpr,
                             logFC)
      })

      results <- data.frame(temp$res,
                            results[setdiff(names(results), c(names(temp$res), names(temp$logres)))],
                            temp$logres)
      results <- reorder.columns(results, start=c("F", "t", "PValue", "P.Value", "FDR"))

      tsmsg("Choosing top DE genes")
      topresults <- results[!is.na(results$FDR) & results$FDR <= 0.1,]
      topresults <- topresults[order(topresults$FDR, topresults$PValue),]
      tsmsgf("Completed limma analysis on %s features in %s samples",
             nrow(feature.counts), ncol(feature.counts))
    })
  }
}

## Same as plotSA, only with ggplot
ggplotSA <- function (fit, xlab = "Average log-expression", ylab = "log2(sigma)",
    zero.weights = FALSE, pch = 16, cex = 0.2, ...)
{
    if (!is(fit, "MArrayLM"))
        stop("fit must be a MArrayLM object")
    x <- fit$Amean
    y <- log2(fit$sigma)
    if (!is.null(fit$weights) && !zero.weights) {
        w <- fit$weights
        w[is.na(w)] <- 0
        w[w < 0] <- 0
        allzero <- apply(w == 0, 1, all)
        y[allzero] <- NA
    }
    (ggplot(data.frame(x,y)) +
     aes(x=x, y=y) +
     geom_point(size=.5) +
     stat_smooth(method="loess", color="red") +
     geom_abline(slope=0, intercept=log2(a$s2.prior)/2, color="blue") +
     ggplot2::xlab(xlab) + ggplot2::ylab(ylab))
}

read.partek <- curry(read.delim, list(sep="\t", na.strings = c("NA", "?", "---")))

do.parteklike.analysis <- function(expdata, feature.counts, gene.lengths=1000, do.logtransform=TRUE, ...) {
  ## Compute RPKM (or RPM)
  tsmsg("Running limma in Partek-like mode.")
  tsmsg("Converting counts to log2(RPKM)")
  libsize <- matrix(nrow=nrow(feature.counts),
                    ncol=ncol(feature.counts),
                    data=colSums(feature.counts),
                    byrow=TRUE)
  glength <- matrix(nrow=nrow(feature.counts),
                    ncol=ncol(feature.counts),
                    data=rep(gene.lengths, length.out=nrow(feature.counts)),
                    byrow=FALSE)
  normfactor <- 1e9 / (glength * libsize)
  feature.RPKM <- normfactor * feature.counts
  if (do.logtransform) {
    feature.RPKM <- log2(1 + feature.RPKM)
  }
  do.limma.analysis(expdata, feature.counts=feature.RPKM, ...,
                    transformation="none", do.ebayes=FALSE)
}

get.fcdesc <- function(fc, log=TRUE) {
  if (log) {
    lfc <- fc
    fc <- 2^lfc
  } else {
    lfc <- log2(fc)
  }
  abslfc <- abs(lfc)
  absfc <- 2^abslfc
  fcsign <- sign(lfc)
  ifelse(fcsign == 0,
         "No change",
         sprintf("%.2f-fold %s", absfc, ifelse(fcsign == 1, "up", "down")))
}

addFPKMToTable <- function(results, gene.lengths) {
  if (!is.null(names(gene.lengths)) && !is.null(rownames(results))) {
    stopifnot(all(rownames(results) %in% names(gene.lengths)))
    gene.lengths <- gene.lengths[rownames(results)]
  }
  stopifnot(length(gene.lengths) == nrow(results))


  cpmcols <- which(str_detect(names(results), "^CPM($|\\.)"))
  suffixes <- str_match(names(results)[cpmcols], "^CPM($|\\..*)")[,2]
  kbmat <- matrix(gene.lengths / 1e3, nrow=nrow(results), ncol=length(suffixes), byrow=FALSE)
  fpkm <- as.matrix(results[cpmcols]) / kbmat
  colnames(fpkm) <- str_c("FPKM", suffixes)

  logcpmcols <- which(str_detect(names(results), "^logCPM($|\\.)"))
  logsuffixes <- str_match(names(results)[logcpmcols], "^logCPM($|\\..*)")[,2]
  logkbmat <- matrix(log2(gene.lengths / 1e3), nrow=nrow(results), ncol=length(logsuffixes), byrow=FALSE)
  logfpkm <- as.matrix(results[logcpmcols]) - logkbmat
  colnames(logfpkm) <- str_c("logFPKM", logsuffixes)

  results <- data.frame(results, results[c(cpmcols, logcpmcols)], check.names=FALSE)
  results[c(cpmcols, logcpmcols)] <- cbind(fpkm, logfpkm)
  names(results)[c(cpmcols, logcpmcols)] <- colnames(cbind(fpkm, logfpkm))
  results
}

drop.zero.columns <- function(x) {
    x[,colSums(abs(x) != 0) > 0]
}
