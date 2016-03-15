source("rnaseq-common.R")
library(DESeq2)

#setMethod("model.matrix", signature=c(object="DataFrame"), definition=function)

get.edger.normfactors <- function(counts, normfactors) {
    if (is.null(normfactors) || is.character(normfactors)) {
        ## A character is interpreted as the method= argument to calcNormFactors
        normfactors <- match.arg(normfactors,
                                 choices=eval(formals(calcNormFactors)$method))
        tsmsgf("Using method %s to calculate sample normalization factors", normfactors)
        normfactors <- calcNormFactors(counts, method=normfactors)
    } else if (is.list(normfactors)) {
        ## A list is interpreted as an arglist for calcNormFactors
        tsmsgf("Using calcNormFactors with to calculate sample normalization factors with the following args")
        ## Print args to stderr compactly
        withOptions(list(max.print=10),
                    Map(message, capture.output(print(normfactors))))
        normfactors$object <- counts
        normfactors <- do.call(calcNormFactors, normfactors)
    }
    else if (is.function(normfactors)) {
        tsmsg("Using user-provided function to calculate sample normalization factors")
        normfactors <- normfactors(counts)
    } else if (is.numeric(normfactors)) {
        tsmsg("Using user-provided sample normalization factors")
    } else {
        stop(str_c("Invalid sample normalization factors:\n", capture.output(str(normfactors, ))))
    }

    if (!is.numeric(normfactors)) {
        stop("Normalization factors nust be numeric")
    } else if (length(normfactors) != ncol(counts)) {
        stop("Need one normalization factor for each column of counts")
    } else if (min(normfactors) <= 0) {
        stop("Need positive normfactors")
    }
    maxlogratio <- abs(log2(Reduce(`/`, range(normfactors))))
    if (maxlogratio >= 1) {
        warning("Normalization factors differ by over 2-fold. Results may be questionable.")
    }
    as.vector(normfactors)
}

get.libsizes <- function(counts, libsizes) {
    ## This works for numeric vectors and NULL, and provides all the
    ## error checking.
    DGEList(counts, lib.size=libsizes)$samples$lib.size
}

check.edger.dimensions <- function(counts, design, feature.annot,
                                   offset, sample.normfactors,
                                   sample.libsizes, feature.lengths,
                                   group, group.blockfactors) {
    if (ncol(counts) <= 1) {
        stop("Need at least 2 columns of counts")
    }

    problems <- character(0)

    ## ncol == nrow
    for (i in c("design", "group.blockfactors")) {
        xi <- get(i)
        tryCatch(stopifnot(is.null(xi) ||
                           ncol(counts) == nrow(xi)),
                 error=function(e) problems <<- c(problems, i))
    }

    ## dim == dim
    for (i in c("offset")) {
        xi <- get(i)
        tryCatch(stopifnot(is.null(xi) ||
                           any(dim(counts) == dim(xi))),
                 error=function(e) problems <<- c(problems, i))
    }

    ## nrow == nrow
    for (i in c("feature.annot")) {
        xi <- get(i)
        tryCatch(stopifnot(is.null(xi) ||
                           nrow(counts) == nrow(xi)),
                 error=function(e) problems <<- c(problems, i))
    }

    ## ncol = length
    for (i in c("sample.normfactors", "sample.libsizes", "group")) {
        xi <- get(i)
        tryCatch(stopifnot(is.null(xi) ||
                           ncol(counts) == length(xi)),
                 error=function(e) problems <<- c(problems, i))
    }

    ## nrow = length
    for (i in c("feature.lengths")) {
        xi <- get(i)
        tryCatch(stopifnot(is.null(xi) ||
                           nrow(counts) == length(xi)),
                 error=function(e) problems <<- c(problems, i))
    }

    if (length(problems) > 0)
        stop(str_c("Aborting analysis due to mismatched dimensions in: ",
                   str_c(problems, collapse=",")))
}

make.group.design <- function(group, block=NULL) {
    group <- droplevels(as.factor(group))
    if (nlevels(group) == 1) {
        groupdesign <- matrix(1, nrow=length(group), ncol=1,
                              dimnames=list(names(group), levels(group)))
    } else {
        groupdesign <- model.matrix(~0+group)
    }
    colnames(groupdesign) <- make.names(levels(group), unique=TRUE, allow_=TRUE)

    if (length(block) == 0)
        return(groupdesign)

    if (is.null(names(block)))
        names(block) <- str_c("block", 1:length(block))
    block <- droplevels(as.data.frame(block))
    ## Use sum-to-zero contrasts for blocking factors, so that the group
    ## coefficients will be the group means after correcting for
    ## blocking factors.
    for (i in names(block)) {
        block[[i]] <- as.factor(block[[i]])
        contrasts(block[[i]]) <- contr.sum(nlevels(block[[i]]))
    }
    blockformula <- as.formula(str_c("~1+", names(block), collapse="+"))
    blockdesign <- model.matrix(blockformula, block)

    attrs <- list(assign=c(attr(groupdesign, "assign"), attr(blockdesign, "assign")[-1] + 1),
                  contrasts=c(attr(groupdesign, "contrasts"), attr(blockdesign, "contrasts")))
    fulldesign <- cbind(groupdesign, blockdesign[,-1])
    attributes(fulldesign) <- c(attrs, attributes(fulldesign))
    colnames(fulldesign) <- make.names(colnames(fulldesign), unique=TRUE, allow_=TRUE)
    fulldesign
}

quantify.by.group.DGEList <- function(dge, group=dge$samples$group, block=NULL, ..., log2=FALSE) {
    stopifnot(!is.null(group))
    group <- droplevels(as.factor(group))
    groupdesign <- make.group.design(group, block)
    groupfit <- glmFit(y=dge, design=groupdesign, ...)
    logCPM <- (groupfit$coefficients[,1:nlevels(group), drop=FALSE] + log(1e+06))/log(2)
    if (log2)
        logCPM
    else
        2 ^ logCPM
}

quantify.by.group.ExpressionSet <- function(eset, group, block=NULL, ..., log2=FALSE) {
    stopifnot(!is.null(group))
    group <- droplevels(as.factor(group))
    groupdesign <- make.group.design(group, block)
    groupfit <- lmFit(object=eset, design=groupdesign, ...)
    logCPM <- groupfit$coefficients[,1:nlevels(group)]/log(2)
    if (log2)
        logCPM
    else
        2 ^ logCPM
}

makeContrast <- function(contrast, levels) {
    makeContrasts(contrasts=list(contrast), levels=levels)[,1]
}

verify.contrast.matrix <- function(contrasts, design) {
    if (!is.numeric(contrasts))
        stop("Contrast matrix must be numeric")
    if (nrow(contrasts) != ncol(design))
        stop("Dimension mismatch between contrasts and design")
    if (ncol(contrasts) >= ncol(design))
        stop("Too many contrasts")
    if (ncol(contrasts) == 0)
        stop("Need at least one contrast")
    if (ncol(contrasts) > 1 && is.null(colnames(contrasts)))
        stop("Multi-contrast spec requires names")
    contrasts
}

make.test.arg <- function(spec, design) {
    if(ncol(design) < 2)
        stop("design must have at least two columns")
    avail.coefs <- colnames(design)
    names(avail.coefs) <- avail.coefs
    contrasts <- NULL
    if (is.matrix(spec)) {
        contrasts <- spec
    } else if (is.character(spec)) {
        ## List of coefs and/or contrasts
        ## Ensure names
        if (is.null(names(spec)))
            names(spec) <- spec
        contrasts <- makeContrasts(contrasts=spec, levels=design)
    } else if (is.numeric(spec)) {
        if (length(spec) > ncol(design)) {
            stop("Contrast has more elements than design has columns")
        }
        else if (length(spec) == ncol(design)) {
            ## Single contrast specified as a numeric vector
            contrasts <- cbind(spec)
        } else {
            ## List of coefs specified as numeric indices
            spec <- avail.coefs[spec]
            contrasts <- makeContrasts(contrasts=spec, levels=design)
        }
    }
    ## edgeR doesn't like a named single-column contrast matrix
    if (ncol(contrasts) == 1)
        colnames(contrasts) <- NULL
    verify.contrast.matrix(contrasts, design)
    contrasts
}

## Basically this is topTags plus adding columns for fold change,
## per-group CPM, and FPKM
get.edger.results.table <- function(dgelrt, groupCPM=NULL, feature.lengths=NULL) {
    table <- data.frame(topTags(dgelrt, n=Inf, sort.by="none"))
    rownames(table) <- rownames(dgelrt)
    table <- table[rownames(dgelrt$coefficients),]

    ## Collect all the rest of the table columns
    table.rest <- table[!str_detect(names(table), "^logFC|logCPM")]

    ## Get log fold changes and log CPM
    logFC <- table[str_detect(names(table), "^logFC")]
    logCPM <- data.frame(table['logCPM'], row.names=rownames(table))
    if (!is.null(groupCPM)) {
        grouplogCPM <- log2(groupCPM)
        groupnames <- str_replace(colnames(groupCPM), "^(CPM\\.?)?", "")
        logCPM <- data.frame(logCPM, grouplogCPM)
        colnames(logCPM) <- str_c("logCPM.", c("mean", groupnames))
    }

    FC <- 2^as.matrix(logFC)
    colnames(FC) <- str_replace(colnames(logFC), "^log", "")
    CPM <- 2^as.matrix(logCPM)
    colnames(CPM) <- str_replace(colnames(logCPM), "^log", "")

    if (!is.null(feature.lengths)) {
        feature.lengths.logkb <- matrix(log2(feature.lengths/1e3), byrow=FALSE,
                                        nrow=nrow(logCPM), ncol=ncol(logCPM))
        logFPKM <- as.matrix(logCPM) - feature.lengths.logkb
        colnames(logFPKM) <- str_replace(colnames(logCPM), "^logCPM", "logFPKM")
        FPKM <- 2^as.matrix(logFPKM)
        colnames(FPKM) <- str_replace(colnames(logFPKM), "^log", "")
    } else {
        logFPKM <- FPKM <- NULL
    }

    allcolumns <- list(row.names=row.names(table),
                       table.rest, logFC, logCPM, logFPKM, FC, CPM, FPKM)
    ## Throw away nulls, because data.frame can't handle them
    allcolumns <- allcolumns[!sapply(allcolumns, is.null)]
    ## Combine into one big data frame
    table <- do.call(data.frame, allcolumns)
    ## Put important columns at the front, unimportant ones at the back
    startcolumns <- c("LR", "F", "PValue", "FDR", colnames(FC))
    endcolumns <- c(colnames(logFC), colnames(logFPKM), colnames(logCPM))
    ## If we have FPKM, then CPM goes at the end
    if (is.null(FPKM)) {
        startcolumns <- c(startcolumns, colnames(CPM))
    } else {
        startcolumns <- c(startcolumns, colnames(FPKM))
        endcolumns <- c(colnames(CPM), endcolumns)
    }
    table <- reorder.columns(table, start=startcolumns, end=endcolumns)
}

topResults <- function(results.table, column="FDR", threshold=0.1, decreasing=FALSE, includeIDs=NULL) {
    full.results.table <- as.data.frame(results.table)
    if (decreasing)
        selected <- !is.na(full.results.table[[column]]) & full.results.table[[column]] >= threshold
    else
        selected <- !is.na(full.results.table[[column]]) & full.results.table[[column]] <= threshold

    results.table <- full.results.table[selected,]
    if (length(includeIDs) > 0) {
        if (is.character(includeIDs))
            extra.results.table <- full.results.table[row.names(full.results.table) %in% includeIDs,]
        else
            extra.results.table <- full.results.table[includeIDs,]
        ## Remove duplicates already in results table
        extra.results.table <- extra.results.table[! row.names(extra.results.table) %in% row.names(results.table),]
        results.table <- rbind(results.table, extra.results.table)
    }
    results.table[order(results.table[[column]], decreasing=decreasing),]
}

do.edger.analysis <-
    function(counts, design,
             tests=NULL,
             feature.annot=NULL,
             feature.lengths=feature.annot$length,
             sample.subset=1:nrow(design),
             feature.subset=1:nrow(counts),
             test.method=c("QL", "LR"),
             offset=NULL,
             sample.normfactors=NULL,
             sample.libsizes=NULL,

             design.disp=design,
             disp.opts=list(),

             group=NULL,
             group.blockfactors=NULL,

             includeIDs=NULL,

             .parallel=TRUE)
{
    stopifnot(nrow(design) == nrow(design.disp))

    test.method <- match.arg(test.method)

    if (!is.null(group.blockfactors))
        group.blockfactors <- as.data.frame(group.blockfactors)

    if (is.null(offset)) {
        sample.normfactors <- get.edger.normfactors(counts, sample.normfactors)
        sample.libsizes <- get.libsizes(counts, sample.libsizes)
    } else {
        tsmsgf("Using user-provided offset matrix.")
        if (!is.null(sample.normfactors) || !is.null(sample.libsizes))
            warning("Ignoring user-provided normalization factors and/or library sizes in favor of offsets.")
        sample.normfactors <- NULL
        sample.libsizes <- NULL
    }

    ## Take appropriate subsets of data and metadata
    if (is.function(sample.subset)) {
        sample.subset <- sample.subset(design)
    }
    if (!is.null(sample.subset)) {
        orig.sample.count <- ncol(counts)
        tsmsg("Subsetting samples")
        counts <- counts[, sample.subset, drop=FALSE]
        design <- design[sample.subset, , drop=FALSE]
        design.disp <- design.disp[sample.subset, , drop=FALSE]
        sample.normfactors <- sample.normfactors[sample.subset]
        sample.libsizes <- sample.libsizes[sample.subset]
        offset <- offset[, sample.subset, drop=FALSE]
        group <- group[sample.subset]
        group.blockfactors <- group.blockfactors[sample.subset, , drop=FALSE]
        tsmsgf("Selected %s out of %s samples.", ncol(counts), orig.sample.count)
    }
    if (is.function(feature.subset)) {
        feature.subset <- feature.subset(counts, feature.annot)
    }
    if (!is.null(feature.subset)) {
        orig.feature.count <- nrow(counts)
        tsmsg("Subsetting features.")
        counts <- counts[feature.subset, , drop=FALSE]
        feature.annot <- feature.annot[feature.subset, , drop=FALSE]
        feature.lengths <- feature.lengths[feature.subset]
        offset <- offset[feature.subset, , drop=FALSE]
        tsmsgf("Selected %s out of %s features.", nrow(counts), orig.feature.count)
    }

    ## Null group is all samples in one group
    if (is.null(group)) {
        group <- rep(1, ncol(counts))
    }
    ## Group must be a factor
    group <- as.factor(group)

    check.edger.dimensions(counts, design, feature.annot, offset, sample.normfactors,
                           sample.libsizes, feature.lengths, group, group.blockfactors)

    tsmsgf("Beginning edgeR analysis on %s features in %s samples",
           nrow(counts), ncol(counts))

    dge <- DGEList(counts, lib.size=sample.libsizes, norm.factors=sample.normfactors,
                   group=group, genes=feature.annot)
    dge$offset <- offset

    tsmsg("Estimating dispersions")
    if (ncol(design.disp) < ncol(design)) {
        tsmsg("Using reduced design matrix for dispersion estimation")
    }
    dge <- do.call(estimateDisp, c(list(y=dge, design=design.disp), disp.opts))

    tsmsg("Fitting GLMs")
    fit <- glmFit(dge, design)
    ## Trended fit for QLFTest
    fit.ql <- glmQLFit(dge, design)

    gof.info <- gof(fit)

    vars.to.return <- c("dge", "fit", "gof.info")

    if (nlevels(group) > 1) {
        tsmsg("Computing group CPM values")
        CPM <- quantify.by.group.DGEList(dge, block=group.blockfactors)
        vars.to.return <- c(vars.to.return, "CPM")
    } else {
        CPM <- NULL
    }

    if (!is.null(tests)) {
        tsmsg("Performing likelihood ratio tests")
        lrts <- llply(tests, function(test) {
            glmLRT(glmfit=fit, contrast=make.test.arg(test, design))
        }, .parallel=.parallel)
        tsmsg("Performing quasi-likelihood F-tests")
        qlfts <- llply(tests, function(test) {
            glmQLFTest(glmfit=fit.ql, contrast=make.test.arg(test, design))
        }, .parallel=.parallel)

        tsmsg("Choosing top DE genes")
        results.lr <- llply(lrts, get.edger.results.table, groupCPM=CPM, feature.lengths=feature.lengths)
        results.ql <- llply(qlfts, get.edger.results.table, groupCPM=CPM, feature.lengths=feature.lengths)
        ## Add GOF outlier column
        results.lr <- llply(results.lr, data.frame, GOF.outlier=gof.info$outlier)
        results.ql <- llply(results.ql, data.frame, GOF.outlier=gof.info$outlier)

        topresults.lr <- llply(results.lr, function(x) {
            x <- topResults(x, includeIDs=includeIDs)
            ordering <- order(x$FDR, x$PValue, -x$LR)
            x[ordering,]
        })
        topresults.ql <- llply(results.ql, function(x) {
            x <- topResults(x, includeIDs=includeIDs)
            ordering <- order(x$FDR, x$PValue, -x$F)
            x[ordering,]
        })

        switch(test.method,
               LR={
                   results <- results.lr
                   topresults <- topresults.lr
               },
               QL={
                   results <- results.ql
                   topresults <- topresults.ql
               },
               stop())

        vars.to.return <- c(vars.to.return, "lrts", "qlfts", "results.lr", "results.ql", "topresults.lr", "topresults.ql", "results", "topresults")
    }

    tsmsgf("Completed %s analysis on %s features in %s samples",
           "edgeR", nrow(counts), ncol(counts))
    mget(vars.to.return, environment())
}

check.deseq.dimensions <- function(counts, design, feature.annot,
                                   feature.lengths, group, group.blockfactors) {
    check.edger.dimensions(counts, design, feature.annot,
                           offset=NULL, sample.normfactors=NULL,
                           sample.libsizes=NULL, feature.lengths,
                           group, group.blockfactors)
}

quantify.by.group.CountDataSet <- function(cds, group=pData(cds)$condition, block=NULL, ..., log2=FALSE) {
    stopifnot(!is.null(group))
    group <- droplevels(as.factor(group))
    groupdesign <- make.group.design(group, block)
    groupfit <- glmFit.CountDataSet(y=cds, design=groupdesign, ...)
    logCPM <- (groupfit$coefficients[,1:nlevels(group)] + log(1e+06))/log(2)
    if (log2)
        logCPM
    else
        2 ^ logCPM
}

do.deseq.analysis <-
    function(counts, design,
             tests=NULL,
             feature.annot=NULL,
             feature.lengths=feature.annot$length,
             sample.subset=1:nrow(design),
             feature.subset=1:nrow(counts),

             disp.args=NULL,
             design.disp=design,

             group=NULL,
             group.blockfactors=NULL,

             includeIDs=NULL,

             .parallel=TRUE)
{
    stopifnot(nrow(design) == nrow(design.disp))

    if (!is.null(group.blockfactors))
        group.blockfactors <- as.data.frame(group.blockfactors)

    ## Take appropriate subsets of data and metadata
    if (is.function(sample.subset)) {
        sample.subset <- sample.subset(design)
    }
    if (!is.null(sample.subset)) {
        orig.sample.count <- ncol(counts)
        tsmsg("Subsetting samples")
        counts <- counts[, sample.subset, drop=FALSE]
        design <- design[sample.subset, , drop=FALSE]
        design.disp <- design.disp[sample.subset, , drop=FALSE]
        group <- group[sample.subset]
        group.blockfactors <- group.blockfactors[sample.subset, , drop=FALSE]
        tsmsgf("Selected %s out of %s samples.", ncol(counts), orig.sample.count)
    }
    if (is.function(feature.subset)) {
        feature.subset <- feature.subset(counts, feature.annot)
    }
    if (!is.null(feature.subset)) {
        orig.feature.count <- nrow(counts)
        tsmsg("Subsetting features.")
        counts <- counts[feature.subset, , drop=FALSE]
        feature.annot <- feature.annot[feature.subset, , drop=FALSE]
        feature.lengths <- feature.lengths[feature.subset]
        tsmsgf("Selected %s out of %s features.", nrow(counts), orig.feature.count)
    }

    ## Null group is all samples in one group
    if (is.null(group)) {
        group <- rep(1, ncol(counts))
    }
    ## Group must be a factor
    group <- as.factor(group)

    check.deseq.dimensions(counts, design, feature.annot,
                           feature.lengths, group, group.blockfactors)

    tsmsgf("Beginning DESeq analysis on %s features in %s samples",
           nrow(counts), ncol(counts))

    cds <- newCountDataSet(counts, conditions=group,
                           featureData=if (!is.null(feature.annot)) as(feature.annot, "AnnotatedDataFrame"))
    tsmsg("Estimating size factors")
    cds <- estimateSizeFactors(cds)

    tsmsg("Estimating dispersions")
    if (ncol(design.disp) < ncol(design)) {
        tsmsg("Using reduced design matrix for dispersion estimation")
    }

    override.disp.args <-
        list(object=cds,
             ## Since we have an already made design instead of a formula
             ## and frame, make a "modelFormula" and "modelFrame" that
             ## simply reproduces the existing design matrix column for
             ## column.
             modelFrame=as.data.frame(design.disp),
             modelFormula=as.formula(str_c("~0+", str_c(colnames(design.disp), collapse="+"))))
    default.disp.args <- list(method="pooled-CR", fitType="parametric",
                              sharingMode="maximum")

    ## Override args take precedence over user args
    disp.args <- merge.defaults.into.list(override.disp.args, disp.args)
    ## Default args only fill in missing values
    disp.args <- merge.defaults.into.list(disp.args, default.disp.args)
    tsmsgf("Using dispersion estimation options:\n%s",
           sprint(unlist(disp.args[c("method", "fitType", "sharingMode")])))

    cds <- do.call(estimateDispersions, disp.args)

    tsmsg("Fitting GLMs")
    fit <- glmFit.CountDataSet(cds, design)
    gof.info <- gof(fit)

    vars.to.return <- c("cds", "fit", "gof.info")

    if (nlevels(group) > 1) {
        tsmsg("Computing group CPM values")
        CPM <- quantify.by.group.CountDataSet(cds, block=group.blockfactors)
        vars.to.return <- c(vars.to.return, "CPM")
    } else {
        CPM <- NULL
    }

    if (!is.null(tests)) {
        tsmsg("Performing likelihood ratio tests")
        lrts <- llply(tests, function(test) {
            glmLRT(glmfit=fit, contrast=make.test.arg(test, design))
        }, .parallel=.parallel)

        tsmsg("Choosing top DE genes")
        results <- llply(lrts, get.edger.results.table, groupCPM=CPM, feature.lengths=feature.lengths)
        ## Add GOF outlier column
        results <- llply(results, data.frame, GOF.outlier=gof.info$outlier)

        topresults <- llply(results, function(x) {
            x <- topResults(x, includeIDs=includeIDs)
            ordering <- order(x$FDR, x$PValue, -x$LR)
            x[ordering,]
        })
        vars.to.return <- c(vars.to.return, "lrts", "results", "topresults")
    }

    tsmsgf("Completed %s analysis on %s features in %s samples",
           "DESeq", nrow(counts), ncol(counts))
    mget(vars.to.return, environment())
}

## Basically this is topTable plus adding columns for fold change,
## per-group CPM, and FPKM
get.limma.results.table <- function(fit, coef, groupCPM=NULL, feature.lengths=NULL) {
    table <- data.frame(topTable(fit, coef=coef, n=Inf, sort.by="none"),
                        row.names=rownames(fit$coefficients))

    ## Change names to be consistent with other analysis pipelines
    renamelist <- c(adj.P.Val="FDR", P.Value="PValue", AveExpr="logCPM")
    ## If there are multiple coefficients, set up renames for their logFC columns
    if (length(coef) > 1) {
        logfc.renames <- setNames(nm=names(coef), sprintf("logFC.%s", names(coef)))
        renamelist <- c(renamelist, logfc.renames)
    } else {
        logfc.renames <- c()
    }
    table <- rename(table, renamelist)

    ## Collect all the rest of the table columns
    table.rest <- table[!str_detect(names(table), "^logFC|logCPM")]

    ## Get log fold changes and log CPM
    logFC <- table[str_detect(names(table), "^logFC")]
    logCPM <- data.frame(table['logCPM'], row.names=rownames(table))
    if (!is.null(groupCPM)) {
        grouplogCPM <- log2(groupCPM)
        groupnames <- str_replace(colnames(groupCPM), "^(CPM\\.?)?", "")
        logCPM <- data.frame(logCPM, grouplogCPM)
        colnames(logCPM) <- str_c("logCPM.", c("mean", groupnames))
    }

    FC <- 2^as.matrix(logFC)
    colnames(FC) <- str_replace(colnames(logFC), "^log", "")
    CPM <- 2^as.matrix(logCPM)
    colnames(CPM) <- str_replace(colnames(logCPM), "^log", "")

    if (!is.null(feature.lengths)) {
        feature.lengths.logkb <- matrix(log2(feature.lengths/1e3), byrow=FALSE,
                                        nrow=nrow(logCPM), ncol=ncol(logCPM))
        logFPKM <- as.matrix(logCPM) - feature.lengths.logkb
        colnames(logFPKM) <- str_replace(colnames(logCPM), "^logCPM", "logFPKM")
        FPKM <- 2^as.matrix(logFPKM)
        colnames(FPKM) <- str_replace(colnames(logFPKM), "^log", "")
    } else {
        logFPKM <- FPKM <- NULL
    }

    allcolumns <- list(row.names=row.names(table),
                       table.rest, logFC, logCPM, logFPKM, FC, CPM, FPKM)
    ## Throw away nulls, because data.frame can't handle them
    allcolumns <- allcolumns[!sapply(allcolumns, is.null)]
    ## Combine into one big data frame
    table <- do.call(data.frame, allcolumns)
    ## Put important columns at the front, unimportant ones at the back
    startcolumns <- c("LR", "F", "PValue", "FDR", colnames(FC))
    endcolumns <- c(colnames(logFC), colnames(logFPKM), colnames(logCPM))
    ## If we have FPKM, then CPM goes at the end
    if (is.null(FPKM)) {
        startcolumns <- c(startcolumns, colnames(CPM))
    } else {
        startcolumns <- c(startcolumns, colnames(FPKM))
        endcolumns <- c(colnames(CPM), endcolumns)
    }
    table <- reorder.columns(table, start=startcolumns, end=endcolumns)
}



voomByGroup <- function(counts, design = NULL, group = NULL, ..., .parallel = FALSE) {
    group <- droplevels(as.factor(group))
    if (nlevels(group) <= 1) {
        return(voom(counts=counts, design=design, ...))
    } else {
        voom.groups <- llply(split(seq_along(group), group), function(lev) {
            voom(counts[, lev, drop = FALSE], design=design[lev, , drop = FALSE], ...)
        }, .parallel=.parallel)
        col.in.main <- seq_along(group)
        col.in.group <- sapply(col.in.main, function(i) sum(group[1:i] == group[i]))
        voom.columns <- llply(col.in.main, function(i) {
            groupi <- group[i]
            coli <- col.in.group[i]
            x <- voom.groups[[groupi]][, coli]
            x$lib.size <- voom.groups[[groupi]]$lib.size[coli]
            x
        })
        ret <- do.call(cbind, unname(voom.columns))
        ret$design <- do.call(rbind, lapply(voom.columns, `[[`, "design"))
        ret$lib.size <- do.call(c, lapply(voom.columns, `[[`, "lib.size"))
        ret$group <- group
        return(ret)
    }
}

transform.by.voom <- function(counts, lib.size = NULL, ...,
                              group = NULL, annot = NULL, .parallel = FALSE) {
    result <- voomByGroup(counts, group=group, lib.size=lib.size, ...)
    if (!is.null(annot)) {
        result$genes <- annot
    }
    result
}
transform.by.vst <- function(counts, lib.size = NULL, ..., annot = NULL) {
    cds <- newCountDataSet(counts, conditions=factor(rep(1, ncol(counts))),
                           featureData=if (!is.null(annot)) as(annot, "AnnotatedDataFrame"))
    if (!is.null(lib.size)) {
        ## Need to convert normalized lib.size to DESeq sizeFactors
        sizeFactors(cds) <- lib.size / geom.mean(colSums(counts(cds)))
    }
    as(cds, "ExpressionSet")
}
transform.by.cpm <- function(counts, lib.size = NULL, norm.factors = NULL, annot = NULL, ...) {
    normdata <- cpm(DGEList(counts, lib.size=lib.size, norm.factors=norm.factors))
    if (is.null(annot)) {
        ExpressionSet(normdata)
    } else {
        ExpressionSet(normdata, featureData=as(annot, "AnnotatedDataFrame"))
    }
}
transform.by.none <- function(counts, annot = NULL, ...) {
    if (is.null(annot)) {
        ExpressionSet(counts)
    } else {
        ExpressionSet(counts, featureData=as(annot, "AnnotatedDataFrame"))
    }
}

do.limma.analysis <-
    function(counts, design,
             tests=NULL,
             feature.annot=NULL,
             feature.lengths=feature.annot$length,
             sample.subset=1:nrow(design),
             feature.subset=1:nrow(counts),
             sample.normfactors=NULL,
             sample.libsizes=NULL,

             transformation=c("voom", "vst", "normcpm", "cpm", "none"),
             voom.group=NULL,
             voom.design=design,

             do.ebayes=TRUE,
             block.design=NULL,

             group=NULL,
             group.blockfactors=NULL,

             includeIDs=NULL,

             .parallel=TRUE)
{
    transformation <- match.arg(transformation)
    if (!is.null(block.design))
        stopifnot(nrow(design) == nrow(block.design))

    if (!is.null(group.blockfactors))
        group.blockfactors <- as.data.frame(group.blockfactors)

    if (transformation == "none") {
        sample.normfactors <- 1
    } else {
        sample.normfactors <- get.edger.normfactors(counts, sample.normfactors)
    }
    sample.libsizes <- get.libsizes(counts, sample.libsizes)

    ## Take appropriate subsets of data and metadata
    if (is.function(sample.subset)) {
        sample.subset <- sample.subset(design)
    }
    if (!is.null(sample.subset)) {
        orig.sample.count <- ncol(counts)
        tsmsg("Subsetting samples")
        counts <- counts[, sample.subset, drop=FALSE]
        design <- design[sample.subset, , drop=FALSE]
        voom.design <- voom.design[sample.subset, , drop=FALSE]
        block.design <- block.design[sample.subset, , drop=FALSE]
        sample.normfactors <- sample.normfactors[sample.subset]
        sample.libsizes <- sample.libsizes[sample.subset]
        group <- group[sample.subset]
        voom.group <- voom.group[sample.subset]
        group.blockfactors <- as.data.frame(group.blockfactors)[sample.subset, , drop=FALSE]
        tsmsgf("Selected %s out of %s samples.", ncol(counts), orig.sample.count)
    }
    if (is.function(feature.subset)) {
        feature.subset <- feature.subset(counts, feature.annot)
    }
    if (!is.null(feature.subset)) {
        orig.feature.count <- nrow(counts)
        tsmsg("Subsetting features.")
        counts <- counts[feature.subset, , drop=FALSE]
        feature.annot <- feature.annot[feature.subset, , drop=FALSE]
        feature.lengths <- feature.lengths[feature.subset]
        tsmsgf("Selected %s out of %s features.", nrow(counts), orig.feature.count)
    }

    ## Null group is all samples in one group
    if (is.null(group)) {
        group <- rep(1, ncol(counts))
    }
    ## Group must be a factor
    group <- as.factor(group)

    norm.lib.sizes <- sample.libsizes * sample.normfactors

    check.edger.dimensions(counts, design, feature.annot, offset=NULL, sample.normfactors,
                           sample.libsizes, feature.lengths, group, group.blockfactors)

    tsmsgf("Beginning limma analysis on %s features in %s samples",
           nrow(counts), ncol(counts))

    eset <-
        switch(transformation,
               voom={
                   tsmsg("Transforming counts to weighted log2(CPM) values using limma voom")
                   transform.by.voom(counts, design=voom.design, group=voom.group,
                                     lib.size=norm.lib.sizes, annot = feature.annot)
               },
               vst={
                   tsmsg("Using DESeq variance-stabilizing transformation")
                   transform.by.vst(counts, lib.size=norm.lib.sizes, annot = feature.annot)
               },
               cpm={
                   tsmsg("Transforming counts to log2(CPM) values")
                   transform.by.cpm(counts, lib.size=norm.lib.sizes, annot = feature.annot)
               },
               none={
                   tsmsg("Skipping transformation step. Assuming input data is already normalized.")
                   transform.by.none(counts, annot = feature.annot)
               })

    if (is.null(block.design)) {
        block.factor <- NULL
        corfit <- NULL
        tsmsg("Fitting linear models")
        fit <- lmFit(eset, design)
    } else {
        tsmsgf("Estimating correlation for blocking factors %s", deparse(unname(block.cols)))
        if (is.factor(block.design)) {
            block.factor <- block.design
        } else {
            block.factor <- designAsFactor(block.design)
        }
        corfit <- duplicateCorrelation(eset, design, block=block.factor)
        tsmsg("Fitting linear models")
        fit <- lmFit(eset, design, block=block.factor, correlation=corfit$consensus.correlation)
    }
    vars.to.return <- c("eset", "fit")

    if (nlevels(group) > 1) {
        tsmsg("Computing group CPM values")
        CPM <- quantify.by.group.ExpressionSet(eset, group=group, block=group.blockfactors, log2=FALSE)
        vars.to.return <- c(vars.to.return, "CPM")
    } else {
        CPM <- NULL
    }

    if (!is.null(tests)) {
        if (do.ebayes) {
            tsmsg("Computing moderated test statistics by eBayes shrinkage")
            eBayes.or.nomod <- eBayes
        } else {
            tsmsg("Computing unmoderated test statistics. This is a bad idea.")
            eBayes.or.nomod <- no.eBayes
        }

        results <- llply(tests, function(test) {
            test <- make.test.arg(test, design)
            if (ncol(test) > 1) {
                if (is.null(colnames(test)))
                    stop("Contrast matrix requires column names")
                coef <- colnames(test)
                names(coef) <- make.names(coef)
            } else {
                coef <- 1
            }
            fit <- contrasts.fit(fit, contrasts=test)
            fit <- eBayes.or.nomod(fit)
            get.limma.results.table(fit, coef, CPM, feature.lengths)
        }, .parallel=.parallel)

        if (transformation == "none")
            warning("Transformation is set to none, so CPM cannot be guaranteed to actually represent counts per million")

        tsmsg("Choosing top DE genes")
        topresults <- llply(results, function(x) {
            topx <- topResults(x, includeIDs=includeIDs)
            topx[order(topx$FDR, topx$PValue),]
        })
        vars.to.return <- c(vars.to.return, "results", "topresults")
    }

    tsmsgf("Completed %s analysis on %s features in %s samples",
           "limma", nrow(counts), ncol(counts))
    mget(vars.to.return, environment())
}

## The following two functions allow the use to edgeR's infrastructure
## to execute the DESeq method. Instead of glmFit(dge, design), use
## glmFit.DESeqSummarizedExperiment(cds, design), then continue as with a normal
## edgeR analysis.
getOffset.DESeqSummarizedExperiment <- function(y) {
  if (any(is.na(sizeFactors(y))))
    stop("Call estimateSizeFactors first")
  log(sizeFactors(y)) - mean(log(sizeFactors(y))) + mean(log(colSums(counts(y))))
}

glmFit.DESeqSummarizedExperiment <- function (y, design = NULL, dispersion = NULL, offset = NULL,
                                 weights = NULL, lib.size = NULL, start = NULL, method = "auto",
                                 ...)
{
  stopifnot(is(y, "DESeqSummarizedExperiment"))
  if (is.null(dispersion)) {
      dispersion <- dispersions(y)
  }
  if (is.null(dispersion)) {
      stop("Call 'estimateDispersions' first before fitting GLMs.")
  }
  if (is.null(offset) && is.null(lib.size))
      offset <- getOffset.DESeqSummarizedExperiment(y)
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
  fit$samples <- as.data.frame(colData(y))
  fit$genes <- as.data.frame(mcols(y))
  new("DGEGLM", fit)
}

quantify.by.group.DESeqSummarizedExperiment <- function(dse, group, block=NULL, ..., log2=FALSE) {
    stopifnot(!is.null(group))
    group <- droplevels(as.factor(group))
    groupdesign <- make.group.design(group, block)
    groupfit <- glmFit.DESeqSummarizedExperiment(y=dse, design=groupdesign, ...)
    logCPM <- (groupfit$coefficients[,1:nlevels(group)] + log(1e+06))/log(2)
    if (log2)
        logCPM
    else
        2 ^ logCPM
}

do.deseq2.analysis <-
    function(counts, design,
             tests=NULL,
             feature.annot=NULL,
             feature.lengths=feature.annot$length,
             sample.subset=1:nrow(design),
             feature.subset=1:nrow(counts),

             disp.fitType="parametric",
             design.disp=design,

             group=NULL,
             group.blockfactors=NULL,

             includeIDs=NULL,

             .parallel=TRUE)
{
    stopifnot(nrow(design) == nrow(design.disp))

    if (!is.null(group.blockfactors))
        group.blockfactors <- as.data.frame(group.blockfactors)

    ## Take appropriate subsets of data and metadata
    if (is.function(sample.subset)) {
        sample.subset <- sample.subset(design)
    }
    if (!is.null(sample.subset)) {
        orig.sample.count <- ncol(counts)
        tsmsg("Subsetting samples")
        counts <- counts[, sample.subset, drop=FALSE]
        design <- design[sample.subset, , drop=FALSE]
        design.disp <- design.disp[sample.subset, , drop=FALSE]
        group <- group[sample.subset]
        group.blockfactors <- group.blockfactors[sample.subset, , drop=FALSE]
        tsmsgf("Selected %s out of %s samples.", ncol(counts), orig.sample.count)
    }
    if (is.function(feature.subset)) {
        feature.subset <- feature.subset(counts, feature.annot)
    }
    if (!is.null(feature.subset)) {
        orig.feature.count <- nrow(counts)
        tsmsg("Subsetting features.")
        counts <- counts[feature.subset, , drop=FALSE]
        feature.annot <- feature.annot[feature.subset, , drop=FALSE]
        feature.lengths <- feature.lengths[feature.subset]
        tsmsgf("Selected %s out of %s features.", nrow(counts), orig.feature.count)
    }

    ## Null group is all samples in one group
    if (is.null(group)) {
        group <- rep(1, ncol(counts))
    }
    ## Group must be a factor
    group <- as.factor(group)

    check.deseq.dimensions(counts, design, feature.annot,
                           feature.lengths, group, group.blockfactors)
    check.deseq.dimensions(counts, design.disp, feature.annot,
                           feature.lengths, group, group.blockfactors)

    tsmsgf("Beginning DESeq2 analysis on %s features in %s samples",
           nrow(counts), ncol(counts))

    ## In DESeq2, "design" refers to the model formula, from which the
    ## design matrix will be constructed. Since we already have a
    ## fully constructed design matrix, we just generate a formula
    ## that simply includes each column of that matrix verbatim.
    model.formula <- as.formula(str_c(c("~0", colnames(design.disp)), collapse="+"))
    dse <- DESeqSummarizedExperimentFromMatrix(countData=counts,
                                               colData=data.frame(design.disp),
                                               design=model.formula)
    tsmsg("Estimating size factors")
    dse <- estimateSizeFactors(dse)

    tsmsg("Estimating dispersions")
    if (ncol(design.disp) < ncol(design)) {
        tsmsg("Using reduced design matrix for dispersion estimation")
    }

    ## This function spams debug messages, so use capture.output to
    ## suppress them unless there's an error
    dsetemp <- NULL
    dse2.output <- capture.output({
        dsetemp <- estimateDispersions(dse, fitType=disp.fitType)
    })
    if (is.null(dsetemp)) {
        message("DESeq2 encountered an error during estimateDispersions. The last few lines of output are shown:")
        cat(tail(dse2.output), sep="\n")
    }
    dse <- dsetemp
    rm(dsetemp)

    tsmsg("Fitting GLMs")
    fit <- glmFit.DESeqSummarizedExperiment(dse, design)
    fit$genes <- feature.annot
    gof.info <- gof(fit)

    vars.to.return <- c("dse", "fit", "gof.info")

    if (nlevels(group) > 1) {
        tsmsg("Computing group CPM values")
        CPM <- quantify.by.group.DESeqSummarizedExperiment(dse, group=group, block=group.blockfactors)
        vars.to.return <- c(vars.to.return, "CPM")
    } else {
        CPM <- NULL
    }

    if (!is.null(tests)) {
        tsmsg("Performing likelihood ratio tests")
        lrts <- llply(tests, function(test) {
            glmLRT(glmfit=fit, contrast=make.test.arg(test, design))
        }, .parallel=.parallel)

        tsmsg("Choosing top DE genes")
        results <- llply(lrts, get.edger.results.table, groupCPM=CPM, feature.lengths=feature.lengths)
        ## Add GOF outlier column
        results <- llply(results, data.frame, GOF.outlier=gof.info$outlier)

        topresults <- llply(results, function(x) {
            x <- topResults(x, includeIDs=includeIDs)
            ordering <- order(x$FDR, x$PValue, -x$LR)
            x[ordering,]
        })
        vars.to.return <- c(vars.to.return, "lrts", "results", "topresults")
    }

    tsmsgf("Completed %s analysis on %s features in %s samples",
           "DESeq2", nrow(counts), ncol(counts))
    mget(vars.to.return, environment())
}
