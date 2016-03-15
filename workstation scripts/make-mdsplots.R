#!/usr/bin/Rscript

library(Biobase)
library(RColorBrewer)
library(Matrix)
library(rtracklayer)
library(GenomicRanges)
library(magrittr)
library(DESeq2)
source("rnaseq-common2.R", chdir=TRUE)

## Having xlxsjars loaded causes a crash
stopifnot(! "package:xlsxjars" %in% search())

llply <- function(...) plyr::llply(..., .paropts=list(.options.multicore=list(preschedule=FALSE, set.seed=FALSE)))

collapse.columns <- function(mat, cols, newname, fun=`+`) {
    ## Convert index to numeric
    cols <- setNames(seq(ncol(mat)), colnames(mat))[cols]
    stopifnot(length(cols) > 0)
    stopifnot(length(cols) < ncol(mat))
    unaffected.cols <- mat[,-cols,drop=FALSE]
    affected.cols <- lapply(cols, function(col) mat[,col])
    new.col <- cbind(Reduce(fun, affected.cols))
    colnames(new.col) <- newname
    cbind(unaffected.cols, new.col)
}

clamp.to.zero <- function(x, tol = 1e-07) {
    ifelse(abs(x) > tol, x, 0)
}

quotemeta <- function(string) {
  str_replace_all(string, "(\\W)", "\\\\\\1")
}

str_detect_prefix <- function(string, prefix) {
    pattern <- perl(str_c("^", quotemeta(prefix)))
    str_detect(string, pattern)
}

str_replace_prefix <- function(string, prefix, ...) {
    pattern <- perl(str_c("^", quotemeta(prefix)))
    str_replace(string=string, pattern=pattern, ...)
}

all.contrasts <- function(f, prefix=lcPrefix(as.character(f))) {
    f <- as.character(unique(f))
    if (length(f) <= 1)
        return(character(0))

    first <- f[1]
    rest <- f[-1]
    mapply(function(x1, x2) {
               ct <- sprintf("(%s) - (%s)", x2, x1)
               if (!is.null(prefix) && str_length(prefix) > 0) {
                   if (all(str_detect_prefix(f, prefix))) {
                       x1 %<>% str_replace_prefix(prefix, "")
                       x2 %<>% str_replace_prefix(prefix, "")
                   }
                   ctname <- sprintf("%s%s.vs.%s",
                                     prefix, x1, x2)
               } else{
                   ctname <- sprintf("%s.vs.%s", x1, x2)
               }
               setNames(ct, make.names(ctname))
           }, first, rest, USE.NAMES=FALSE)
}

all.pairwise.contrasts <- function(f, prefix=lcPrefix(as.character(f))) {
    f <- as.character(unique(f))
    if (length(f) <= 1)
        return(character(0))

    f %>% combn(2, function(x) {
                    x1 <- x[1]
                    x2 <- x[2]
                    ct <- sprintf("(%s) - (%s)", x2, x1)
                    if (!is.null(prefix) && str_length(prefix) > 0) {
                        if (all(str_detect_prefix(f, prefix))) {
                            x1 %<>% str_replace_prefix(prefix, "")
                            x2 %<>% str_replace_prefix(prefix, "")
                        }
                        ctname <- sprintf("%s%s.vs.%s",
                                          prefix, x1, x2)
                    } else{
                        ctname <- sprintf("%s.vs.%s", x1, x2)
                    }
                    setNames(ct, make.names(ctname))
                }, simplify = FALSE) %>% unlist
}

coef.mean <- function(coef.names) {
    coef.names <- unique(coef.names)
    coefs.sum <- str_join(coef.names, collapse=" + ")
    sprintf("((%s) / %s)", coefs.sum, length(coef.names))
}

drop.constant.columns <- function(mat, intcol=1, tol=1e-6) {
    if (min(dim(mat)) == 0) {
        mat
    } else {
        intcol <- setNames(seq(ncol(mat)), colnames(mat))[intcol]
        colRanges <- rowMax(t(mat)) - rowMin(t(mat))
        remcols <- !seq(ncol(mat)) %in% intcol & colRanges < tol
        mat[,!remcols, drop=FALSE]
    }
}

subdesign <- function(design, subset, intcol=str_detect(colnames(design), "Intercept")) {
    subdes <- design[subset, , drop=FALSE]
    subdes <- drop.constant.columns(subdes, intcol=intcol)
    if (rankMatrix(subdes) < ncol(subdes)) {
        warning("Subdesign matrix is not full rank")
    }
    if (rankMatrix(subdes) >= nrow(subdes)) {
        warning("Subdesign matrix is has no residual df")
    }
    subdes
}

addsummary <- function(tables) {
    summary <- data.frame(diffcount=sapply(tables, nrow))
    c(list(summary=summary), tables)
}

flatten.CharacterList <- function(x, sep=",") {
  if (is.list(x)) {
    x[!is.na(x)] <- laply(x[!is.na(x)], str_c, collapse=sep, .parallel=TRUE)
    x <- as(x, "character")
  }
  x
}

fix.cl.in.df <- function(df) {
  data.frame(llply(df, flatten.CharacterList, .parallel=TRUE), row.names=row.names(df))
}

group.by.factor <- function(f) split(seq_along(f), f)

suppressPlot <- function(arg) {
    dev.new()
    result <- arg
    dev.off()
    result
}

suppressPlot <- function(arg) {
    png("/dev/null")
    on.exit(dev.off())
    result <- arg
    result
}

add.numbered.colnames <- function(x, prefix="C") {
    x %>% set_colnames(sprintf("%s%i", prefix, seq(from=1, length.out=ncol(x))))
}

## Version of cpm that uses an offset matrix instead of lib sizes
cpmWithOffset <- function(dge, offset=expandAsMatrix(getOffset(dge), dim(dge)),
                          log = FALSE, prior.count = 0.25, ...) {
    x <- dge$counts
    effective.lib.size <- exp(offset)
    if (log) {
        prior.count.scaled <- effective.lib.size/mean(effective.lib.size) * prior.count
        effective.lib.size <- effective.lib.size + 2 * prior.count.scaled
    }
    effective.lib.size <- 1e-06 * effective.lib.size
    if (log)
        log2((x + prior.count.scaled) / effective.lib.size)
    else x / effective.lib.size
}

## Like limma's removeBatchEffect, but for DGEList. Modifies the
## offsets instead of the data.
offsetBatchEffect <- function (dge, batch = NULL, batch2 = NULL, covariates = NULL,
                               design = matrix(1, ncol(dge), 1), ...) {
    if (is.null(batch) && is.null(batch2) && is.null(covariates))
        return(dge)
    if (!is.null(batch)) {
        batch <- as.factor(batch)
        contrasts(batch) <- contr.sum(levels(batch))
        batch <- model.matrix(~batch)[, -1, drop = FALSE]
    }
    if (!is.null(batch2)) {
        batch2 <- as.factor(batch2)
        contrasts(batch2) <- contr.sum(levels(batch2))
        batch2 <- model.matrix(~batch2)[, -1, drop = FALSE]
    }
    if (!is.null(covariates))
        covariates <- as.matrix(covariates)
    X.batch <- cbind(batch, batch2, covariates)
    fit <- glmFit(dge, cbind(design, X.batch), ..., priot.count = 0)
    beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
    beta[is.na(beta)] <- 0
    dge$ofsfset <- fit$offset + beta %*% t(X.batch)
    dge
}

## Version of voom that uses an offset matrix instead of lib sizes
voomWithOffset <- function (dge, design = NULL, offset=expandAsMatrix(getOffset(dge), dim(dge)),
                            normalize.method = "none", plot = FALSE, span = 0.5, ...)
{
    out <- list()
    out$genes <- dge$genes
    out$targets <- dge$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) >
        0)
        design <- model.matrix(~group, data = counts$samples)
    counts <- dge$counts
    if (is.null(design)) {
        design <- matrix(1, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
    }

    effective.lib.size <- exp(offset)

    y <- log2((counts + 0.5)/(effective.lib.size + 1) * 1e+06)
    y <- normalizeBetweenArrays(y, method = normalize.method)
    fit <- lmFit(y, design, ...)
    if (is.null(fit$Amean))
        fit$Amean <- rowMeans(y, na.rm = TRUE)
    sx <- fit$Amean + mean(log2(effective.lib.size + 1)) - log2(1e+06)
    sy <- sqrt(fit$sigma)
    allzero <- rowSums(counts) == 0
    if (any(allzero)) {
        sx <- sx[!allzero]
        sy <- sy[!allzero]
    }
    l <- lowess(sx, sy, f = span)
    if (plot) {
        plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )",
            pch = 16, cex = 0.25)
        title("voom: Mean-variance trend")
        lines(l, col = "red")
    }
    f <- approxfun(l, rule = 2)
    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,
            j, drop = FALSE])
    }
    else {
        fitted.values <- fit$coef %*% t(fit$design)
    }
    fitted.cpm <- 2^fitted.values
    ## fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.count <- 1e-06 * fitted.cpm * (effective.lib.size + 1)
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    out$E <- y
    out$weights <- w
    out$design <- design
    out$effective.lib.size <- effective.lib.size
    if (is.null(out$targets))
        out$targets <- data.frame(lib.size = exp(colMeans(offset)))
    else out$targets$lib.size <- exp(colMeans(offset))
    new("EList", out)
}


## Version of voom that uses an offset matrix instead of lib sizes
voomWithQualityWeightsAndOffset <-
    function (dge, design = NULL,
              offset=expandAsMatrix(getOffset(dge), dim(dge)),
              normalize.method = "none",
              plot = FALSE, span = 0.5, var.design = NULL, method = "genebygene",
              maxiter = 50, tol = 1e-10, trace = FALSE, replace.weights = TRUE,
              col = NULL, ...)
{
    counts <- dge$counts
    if (plot) {
        oldpar <- par(mfrow = c(1, 2))
        on.exit(par(oldpar))
    }
    v <- voomWithOffset(dge, design = design, offset = offset, normalize.method = normalize.method,
        plot = FALSE, span = span, ...)
    aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, var.design = var.design)
    v <- voomWithOffset(dge, design = design, weights = aw, offset = offset,
        normalize.method = normalize.method, plot = plot, span = span,
        ...)
    aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, trace = trace, var.design = var.design)
    wts <- asMatrixWeights(aw, dim(v)) * v$weights
    attr(wts, "arrayweights") <- NULL
    if (plot) {
        barplot(aw, names = 1:length(aw), main = "Sample-specific weights",
            ylab = "Weight", xlab = "Sample", col = col)
        abline(h = 1, col = 2, lty = 2)
    }
    if (replace.weights) {
        v$weights <- wts
        v$sample.weights <- aw
        return(v)
    }
    else {
        return(wts)
    }
}

estimateDispByGroup <- function(dge, group=as.factor(dge$samples$group), batch, ...) {
    stopifnot(nlevels(group) > 1)
    stopifnot(length(group) == ncol(dge))
    if (!is.list(batch)) {
        batch <- list(batch=batch)
    }
    batch <- as.data.frame(batch)
    stopifnot(nrow(batch) == ncol(dge))
    colnames(batch) %>% make.names(unique=TRUE)
    igroup <- seq_len(ncol(dge)) %>% split(group)
    lapply(igroup, function(i) {
        group.dge <- dge[,i]
        group.batch <- droplevels(batch[i,, drop=FALSE])
        group.batch <- group.batch[sapply(group.batch, . %>% unique %>% length %>% is_greater_than(1))]
        group.vars <- names(group.batch)
        if (length(group.vars) == 0)
            group.vars <- "1"
        group.batch.formula <- as.formula(str_c("~", str_c(group.vars, collapse="+")))
        des <- model.matrix(group.batch.formula, group.batch)
        estimateDisp(group.dge, des, ...)
    })
}

sexp <- readRDS("promoter-group-counts-1kb.RDS") %>% updateObject

## expdata <- droplevels(as.data.frame(colData(sexp)))
expdata <- read.xlsx("sampledata.xlsx", 1)
rownames(expdata) <- expdata$Sample
expdata$Lane <- factor(expdata$Lane)
expdata$Timepoint <- factor(expdata$Timepoint, levels=str_c("D", c(0, 1, 5, 14)))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
expdata$Sampletype <- factor(expdata$Sampletype, levels=c("input", "H3K4me2", "H3K4me3", "H3K27me3"))
## Ensure expdata has same sample order as sexp
expdata <- expdata[colData(sexp)$Sample,]
featuredata <- as.data.frame(mcols(sexp))
featuredata <- fix.cl.in.df(featuredata)
stopifnot(!any(sapply(featuredata, is.list)))
counts <- assays(sexp)$counts

## Combine technical replicates
expdata$Biolrep <- with(expdata, interaction(Donor, Celltype, Sampletype, Timepoint, sep=".", drop=TRUE))
counts <- do.call(cbind, llply(split(seq(nrow(expdata)), expdata$Biolrep), function(i) {
    setNames(rowSums(counts[,i,drop=FALSE]),
             rownames(counts))
}, .parallel=TRUE))
expdata <- expdata[!duplicated(expdata$Biolrep),
                   setdiff(colnames(expdata), c("Sample", "Run", "Lane", "BAM"))]
rownames(expdata) <- expdata$Biolrep
expdata <- expdata[colnames(counts),]

expdata %<>%
    mutate(
        Condition=interaction(Celltype, Timepoint, sep="."),
        Group=interaction(Sampletype, Condition, sep=".")
    )

outliers <- c("4659.Naive.H3K4me3.D1",
              "5291.Memory.H3K27me3.D14",
              "5291.Naive.input.D1",
              "5291.Memory.input.D1",
              "4659.Naive.input.D1",
              "4659.Memory.input.D1",
              "5291.Naive.input.D0")

condcolors <- suppressWarnings(brewer.pal(Inf, "Paired")) %>%
    matrix(nrow=2) %>% .[,c(3,2,1,5)] %>% as.vector %>%
    setNames(levels(expdata$Condition))

## Solid
donorshapes <- c(16, 17, 15, 18)
## Outlines
## donorshapes <- c(1, 2, 0, 5)

{
    tsmsg("Doing standard edgeR PCoA plots")
    mds.points <- lapply(c("input", "H3K4me2", "H3K4me3") %>% setNames(nm=.), function(i) {
        icols <- expdata$Sampletype == i & ! (expdata$Biolrep  %in% outliers)
        subcounts <- counts[,icols]
        subexpdata <- expdata[icols,] %>% droplevels
        dge <- DGEList(subcounts) %>% calcNormFactors
        dmat <- suppressPlot(dge %>% plotMDS %$% distance.matrix %>% as.dist)
        mds <- cmdscale(dmat, k=attr(dmat, "Size") - 1, eig=TRUE)
        mds$points %>% add.numbered.colnames("PC") %>% data.frame(subexpdata, .)
    })
    mds.group.mean.points <- lapply(mds.points, function(x) {
        x %>%
            group_by(Group) %>%
            do({
                df <- .
                dimcols <- str_detect(colnames(df), "^PC\\d+$")
                means <- rbind(colMeans(as.matrix(.[,dimcols])))
                df[1,!dimcols] %>% cbind(means)
            })
    })
    local({
        pdf("results/1kb/mdsplots-pubready.pdf")
        on.exit(dev.off())
        for (i in names(mds.points)) {
            print(ggplot(mds.points[[i]]) +
                  aes(x=PC1, y=PC2, color=Condition, fill=Condition, shape=Donor, group=Celltype) +
                  geom_point(size=3) +
                  geom_point(shape=10, size=5, data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  geom_path(data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  coord_equal() +
                  scale_color_manual(values=condcolors) +
                  scale_fill_manual(values=condcolors) +
                  scale_shape_manual(values=donorshapes) +
                  ggtitle(sprintf("PCoA plot for %s samples (PCs 1 & 2)", i)))
            print(ggplot(mds.points[[i]]) +
                  aes(x=PC3, y=PC4, color=Condition, fill=Condition, shape=Donor, group=Celltype) +
                  geom_point(size=3) +
                  geom_point(shape=10, size=5, data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  geom_path(data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  coord_equal() +
                  scale_color_manual(values=condcolors) +
                  scale_fill_manual(values=condcolors) +
                  scale_shape_manual(values=donorshapes) +
                  ggtitle(sprintf("PCoA plot for %s samples (PCs 3 & 4)", i)))
        }
    })
}

{
    tsmsg("Doing edgeR PCoA with donor batch correction")
    mds.points <- lapply(c("input", "H3K4me2", "H3K4me3") %>% setNames(nm=.), function(i) {
        icols <- expdata$Sampletype == i & ! (expdata$Biolrep  %in% outliers)
        subcounts <- counts[,icols]
        subexpdata <- expdata[icols,] %>% droplevels
        dge <- DGEList(subcounts) %>% calcNormFactors
        design <- model.matrix(~Group + Donor, subexpdata)
        dge %<>%
            estimateDisp(design, robust=TRUE) %>%
            offsetBatchEffect(design=model.matrix(~Group, subexpdata), batch=subexpdata$Donor)
        dmat <- suppressPlot(
            dge %>% cpmWithOffset(log=TRUE, prior.count=2) %>%
                plotMDS %$% distance.matrix %>% as.dist)
        mds <- cmdscale(dmat, k=attr(dmat, "Size") - 1, eig=TRUE)
        mds$points %>% add.numbered.colnames("PC") %>% data.frame(subexpdata, .)
    })
    mds.group.mean.points <- lapply(mds.points, function(x) {
        x %>%
            group_by(Group) %>%
            do({
                df <- .
                dimcols <- str_detect(colnames(df), "^PC\\d+$")
                means <- rbind(colMeans(as.matrix(.[,dimcols])))
                df[1,!dimcols] %>% cbind(means)
            })
    })
    local({
        pdf("results/1kb/mdsplots-pubready-bc.pdf")
        on.exit(dev.off())
        for (i in names(mds.points)) {
            print(ggplot(mds.points[[i]]) +
                  aes(x=PC1, y=PC2, color=Condition, fill=Condition, shape=Donor, group=Celltype) +
                  geom_point(size=3) +
                  geom_point(shape=10, size=5, data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  geom_path(data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  coord_equal() +
                  scale_color_manual(values=condcolors) +
                  scale_fill_manual(values=condcolors) +
                  scale_shape_manual(values=donorshapes) +
                  ggtitle(sprintf("PCoA plot for %s samples (PCs 1 & 2)", i)))
            print(ggplot(mds.points[[i]]) +
                  aes(x=PC3, y=PC4, color=Condition, fill=Condition, shape=Donor, group=Celltype) +
                  geom_point(size=3) +
                  geom_point(shape=10, size=5, data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  geom_path(data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  coord_equal() +
                  scale_color_manual(values=condcolors) +
                  scale_fill_manual(values=condcolors) +
                  scale_shape_manual(values=donorshapes) +
                  ggtitle(sprintf("PCoA plot for %s samples (PCs 3 & 4)", i)))
        }
    })
}

{
    tsmsg("Doing standard voom PCoA")
    mds.points <- lapply(c("input", "H3K4me2", "H3K4me3") %>% setNames(nm=.), function(i) {
        icols <- expdata$Sampletype == i & ! (expdata$Biolrep  %in% outliers)
        subcounts <- counts[,icols]
        subexpdata <- expdata[icols,] %>% droplevels
        dge <- DGEList(subcounts) %>% calcNormFactors
        design <- model.matrix(~Group + Donor, subexpdata)
        elist <- voom(dge, design)
        dmat <- suppressPlot(elist %>% plotMDS %$% distance.matrix %>% as.dist)
        mds <- cmdscale(dmat, k=attr(dmat, "Size") - 1, eig=TRUE)
        mds$points %>% add.numbered.colnames("PC") %>% data.frame(subexpdata, .)
    })
    mds.group.mean.points <- lapply(mds.points, function(x) {
        x %>%
            group_by(Group) %>%
            do({
                df <- .
                dimcols <- str_detect(colnames(df), "^PC\\d+$")
                means <- rbind(colMeans(as.matrix(.[,dimcols])))
                df[1,!dimcols] %>% cbind(means)
            })
    })
    local({
        pdf("results/1kb/mdsplots-pubready-voom.pdf")
        on.exit(dev.off())
        for (i in names(mds.points)) {
            print(ggplot(mds.points[[i]]) +
                  aes(x=PC1, y=PC2, color=Condition, fill=Condition, shape=Donor, group=Celltype) +
                  geom_point(size=3) +
                  geom_point(shape=10, size=5, data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  geom_path(data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  coord_equal() +
                  scale_color_manual(values=condcolors) +
                  scale_fill_manual(values=condcolors) +
                  scale_shape_manual(values=donorshapes) +
                  ggtitle(sprintf("PCoA plot for %s samples (PCs 1 & 2)", i)))
            print(ggplot(mds.points[[i]]) +
                  aes(x=PC3, y=PC4, color=Condition, fill=Condition, shape=Donor, group=Celltype) +
                  geom_point(size=3) +
                  geom_point(shape=10, size=5, data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  geom_path(data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  coord_equal() +
                  scale_color_manual(values=condcolors) +
                  scale_fill_manual(values=condcolors) +
                  scale_shape_manual(values=donorshapes) +
                  ggtitle(sprintf("PCoA plot for %s samples (PCs 3 & 4)", i)))
        }
    })
}

{
    tsmsg("Doing voom PCoA with donor batch correction")
    mds.points <- lapply(c("input", "H3K4me2", "H3K4me3") %>% setNames(nm=.), function(i) {
        icols <- expdata$Sampletype == i & ! (expdata$Biolrep  %in% outliers)
        subcounts <- counts[,icols]
        subexpdata <- expdata[icols,] %>% droplevels
        dge <- DGEList(subcounts) %>% calcNormFactors
        design <- model.matrix(~Group + Donor, subexpdata)
        dge %<>% estimateDisp(design, robust=TRUE) %>%
            offsetBatchEffect(design=model.matrix(~Group, subexpdata), batch=subexpdata$Donor)
        elist <- voomWithOffset(dge, design)
        dmat <- suppressPlot(elist %>% plotMDS %$% distance.matrix %>% as.dist)
        mds <- cmdscale(dmat, k=attr(dmat, "Size") - 1, eig=TRUE)
        mds$points %>% add.numbered.colnames("PC") %>% data.frame(subexpdata, .)
    })
    mds.group.mean.points <- lapply(mds.points, function(x) {
        x %>%
            group_by(Group) %>%
            do({
                df <- .
                dimcols <- str_detect(colnames(df), "^PC\\d+$")
                means <- rbind(colMeans(as.matrix(.[,dimcols])))
                df[1,!dimcols] %>% cbind(means)
            })
    })
    local({
        pdf("results/1kb/mdsplots-pubready-voom-bc.pdf")
        on.exit(dev.off())
        for (i in names(mds.points)) {
            print(ggplot(mds.points[[i]]) +
                  aes(x=PC1, y=PC2, color=Condition, fill=Condition, shape=Donor, group=Celltype) +
                  geom_point(size=3) +
                  geom_point(shape=10, size=5, data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  geom_path(data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  coord_equal() +
                  scale_color_manual(values=condcolors) +
                  scale_fill_manual(values=condcolors) +
                  scale_shape_manual(values=donorshapes) +
                  ggtitle(sprintf("PCoA plot for %s samples (PCs 1 & 2)", i)))
            print(ggplot(mds.points[[i]]) +
                  aes(x=PC3, y=PC4, color=Condition, fill=Condition, shape=Donor, group=Celltype) +
                  geom_point(size=3) +
                  geom_point(shape=10, size=5, data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  geom_path(data=mds.group.mean.points[[i]], show_guide=FALSE) +
                  coord_equal() +
                  scale_color_manual(values=condcolors) +
                  scale_fill_manual(values=condcolors) +
                  scale_shape_manual(values=donorshapes) +
                  ggtitle(sprintf("PCoA plot for %s samples (PCs 3 & 4)", i)))
        }
    })
}
