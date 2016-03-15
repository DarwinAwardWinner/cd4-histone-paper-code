#!/usr/bin/Rscript

library(RColorBrewer)
library(Matrix)
source("rnaseq-common2.R", chdir=TRUE)
library(statmod)
library(GMD)
library(amap)

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

all.contrasts <- function(f) {
    f <- as.character(unique(f))
    if (length(f) <= 1)
        return(character(0))
    else
        sprintf("(%s) - (%s)", f[-1], f[1])
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

suppressPlot <- function(expr) {
    pdf("/dev/null")
    on.exit(dev.off())
    expr
}

mds.dist <- function(...) {
    as.dist(suppressPlot(plotMDS(...)$distance.matrix))
}

brewer.pal.full <- function(name) {
    suppressWarnings(brewer.pal(Inf, name))
}

dir.create("results", showWarnings=FALSE)

bigbin.sexp <- readRDS("bigbin-counts.RDS")
sexps <- readRDS("peak-counts.RDS")

## Throw away ChIP input samples, calculate norm factors from bigbins
sexps <- llply(sexps, function(x) {
    x <- x[,colData(x)$Sampletype != "input"]
    bbx <- bigbin.sexp[,colnames(x)]
    colData(x)$NormFactor <- calcNormFactors(assays(bbx)$counts)
    colnames(x) <- with(colData(x), interaction(Celltype, Timepoint, Donor))
    x
})

known.outliers <- list(H3K27me3=c("Memory.D14.5291"))

heatmap.colscheme <- colorRampPalette(rev(brewer.pal.full("Greens")))

edger.analyses <- llply(names(sexps), function(chip) {
    sexp <- sexps[[chip]]
    sexp <- sexp[,!colnames(sexp) %in% known.outliers[[chip]]]
    expdata <- as(colData(sexp), "data.frame")
    expdata$Timepoint <- factor(expdata$Timepoint, levels=str_c("D", c(0, 1, 5, 14)))
    expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
    design <- make.group.design(group=with(expdata, Celltype:Timepoint), block=expdata$Donor)
    all.timepoints <- levels(expdata$Timepoint)
    nonzero.timepoints <- setdiff(all.timepoints, "D0")
    timepoint.anova.tests <- foreach(ctype=c("Naive", "Memory"), .combine=c) %do% {
        x <- sprintf("%s.%s - %s.D0", ctype, nonzero.timepoints, ctype)
        names(x) <- sprintf("%s.D0v%s", ctype, nonzero.timepoints)
        setNames(nm=sprintf("%s.AllT", ctype), list(x))
    }
    timepoint.individual.tests <- unlist(unname(timepoint.anova.tests))

    celltype.singlet.tests <- sprintf("Memory.%s - Naive.%s",
                                      all.timepoints, all.timepoints)
    names(celltype.singlet.tests) <- sprintf("NvM.%s", all.timepoints)
    celltype.allt.test <- list(NvM.AllT=celltype.singlet.tests)

    factorial.singlet.tests <- sprintf("(Memory.%s - Naive.%s) - (Memory.D0 - Naive.D0)", nonzero.timepoints, nonzero.timepoints)
    names(factorial.singlet.tests) <- sprintf("Fac.%s", nonzero.timepoints)
    factorial.allt.test <- list(Fac.AllT=factorial.singlet.tests)

    mi.vs.nf.test <- c(ND14vMD0="Memory.D0 - Naive.D14")

    alltests <- c(timepoint.anova.tests,
                  timepoint.individual.tests,
                  celltype.allt.test,
                  celltype.singlet.tests,
                  factorial.allt.test,
                  factorial.singlet.tests,
                  mi.vs.nf.test)
    alltest.contrasts <- lapply(alltests, function(test) makeContrasts(contrasts=test, levels=design))
    dge <- DGEList(assays(sexp)$counts, group=with(expdata, Celltype:Timepoint), norm.factors=expdata$NormFactor)
    ## dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design)
    qlfts <- llply(alltest.contrasts, function(cont) glmQLFTest(dge, design, contrast=cont, robust=TRUE), .parallel=TRUE)
    tts <- llply(qlfts, function(x) topTags(x, n=Inf, sort.by="PValue")$table)
    ## All significant at FDR=0.2 or top 100
    toppeaks <- llply(tts, function(tt) tt[tt$FDR <= 0.2 | seq(nrow(tt)) <= 50,])
    subdges <- c(list(All=dge),
                 llply(toppeaks, function(x) dge[rownames(x),]))
    dmats.mds <- llply(subdges, mds.dist)
    dmats.cor <- llply(subdges, function(dge) gdist(t(cpm(dge, log=TRUE)), method="correlation.of.observations"))
    dmats.pear <- llply(subdges, function(dge) Dist(t(cpm(dge, log=TRUE)), method="pearson"))
    pdf(sprintf("results/Peak MDSDist Heatmaps %s.pdf", chip))
    for (i in names(dmats.mds)) {
        heatmap.3(dmats.mds[[i]], color.FUN=heatmap.colscheme, breaks=256,
                  dendrogram="both", Colv=TRUE, Rowv=TRUE, revC=TRUE,
                  main=sprintf("%s Peak Heatmap (MDS distances, %s peaks)", chip, i),
                  cex.main=1)
    }
    dev.off()
    pdf(sprintf("results/Peak CorDist Heatmaps %s.pdf", chip))
    for (i in names(dmats.cor)) {
        heatmap.3(dmats.cor[[i]], color.FUN=heatmap.colscheme, breaks=256,
                  dendrogram="both", Colv=TRUE, Rowv=TRUE, revC=TRUE,
                  main=sprintf("%s Peak Heatmap (Correlation distances, %s peaks)", chip, i),
                  cex.main=1)
    }
    dev.off()
    pdf(sprintf("results/Peak PearsonDist Heatmaps %s.pdf", chip))
    for (i in names(dmats.pear)) {
        heatmap.3(dmats.pear[[i]], color.FUN=heatmap.colscheme, breaks=256,
                  dendrogram="both", Colv=TRUE, Rowv=TRUE, revC=TRUE,
                  main=sprintf("%s Peak Heatmap (Pearson distances, %s peaks)", chip, i),
                  cex.main=1)
    }
    dev.off()
    list(dge=dge,
         design=design,
         alltests=alltests,
         qlfts=qlfts,
         tts=llply(qlfts, function(x) topTags(x, n=Inf, sort.by="none")$table))
})
