#!usr/bin/env Rscript

library(RColorBrewer)
library(Matrix)
library(rtracklayer)
library(GenomicRanges)
source("rnaseq-common2.R", chdir=TRUE)
library(magrittr)
library(dplyr)
library(annotate)
library(org.Hs.eg.db)

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

read.and.collapse.techreps <- function(filename, expdata) {
    sexp <- readRDS(filename)

    ## ## expdata <- droplevels(as.data.frame(colData(sexp)))
    ## expdata <- read.xlsx("sampledata.xlsx", 1)
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

    expdata$Group <- with(expdata, {
        x <- as.character(Sampletype:Celltype:Timepoint)
        ## ## Combine all input samples into one group
        ## x[Sampletype == "input"] <- "input"
        factor(x)
    })

    newsexp <- SummarizedExperiment(
        assays=list(counts=counts),
        rowData=rowData(sexp),
        colData=DataFrame(expdata))
    mcols(newsexp) <- featuredata
    newsexp
}

fixup.mcols <- function(sexp) {
    entrez <- sexp %>% mcols %$% ENTREZ %>% as.character
    fixed.entrez.annot <- c("SYMBOL", "GENENAME", "UNIGENE",
                            "ENSEMBL", "REFSEQ", "UNIPROT") %>%
        setNames(llply(., . %>% {
                  CharacterList(lookUp(entrez,
                                       data="org.Hs.eg",
                                       what=.))
              }), .)
    mcols(sexp)[names(fixed.entrez.annot)] <- fixed.entrez.annot[]
    sexp
}

split.and.compute.sample.FPKMs <- function(sexp, group) {
    splits <- split(seq_along(group), group, drop=TRUE)
    lapply(splits, function(i) {
        subsexp <- sexp[,i]
        cts <- subsexp %>% assay("counts")
        cpms <- cts %>% DGEList %>% calcNormFactors(method="TMM") %>%
            cpm()
        promoter.widths.in.kb <- subsexp %>% rowData %>% reduce %>%
            width %>% sum %>% divide_by(1000)
        fpkms <- cpms / promoter.widths.in.kb
        assay(subsexp, "CPM") <- cpms
        assay(subsexp, "FPKM") <- fpkms
        subsexp
    })
}

expdata <- read.xlsx("sampledata.xlsx", 1)

sexp1 <- read.and.collapse.techreps("promoter-group-counts-1kb.RDS", expdata)
sexp1 %<>% fixup.mcols
sexp2.5 <- read.and.collapse.techreps("promoter-group-counts-2.5kb.RDS", expdata)
sexp2.5 %<>% fixup.mcols

s1.by.stype <- split.and.compute.sample.FPKMs(sexp1, colData(sexp1)$Sampletype)
s1.sampleFPKM.by.stype <- lapply(s1.by.stype, . %>% assay("FPKM"))
s1.groupFPKM.by.stype <- lapply(s1.by.stype, . %>% {
    cts <- assay(., "counts")
    expdata <- {.} %>% colData %>% as("data.frame") %>% droplevels
    design <- make.group.design(expdata$Group, expdata["Donor"])
    dge <- cts %>% DGEList %>% calcNormFactors %>% estimateDisp(design=design, robust=TRUE)
    quantify.by.group.DGEList(dge, expdata$Group, expdata["Donor"])
})
s1.sampleFPKM.by.stype %<>% lapply(. %>% data.frame %>% {
    cbind(mcols(sexp1)[c("ENTREZ", "SYMBOL", "GENENAME")] %>% as("data.frame"), .)
})
s1.groupFPKM.by.stype %<>% lapply(. %>% data.frame %>% {
    cbind(mcols(sexp1)[c("ENTREZ", "SYMBOL", "GENENAME")] %>% as("data.frame"), .)
})

s2.5.by.stype <- split.and.compute.sample.FPKMs(sexp2.5, colData(sexp2.5)$Sampletype)
s2.5.sampleFPKM.by.stype <- lapply(s2.5.by.stype, . %>% assay("FPKM"))
s2.5.groupFPKM.by.stype <- lapply(s2.5.by.stype, . %>% {
    cts <- assay(., "counts")
    expdata <- {.} %>% colData %>% as("data.frame") %>% droplevels
    design <- make.group.design(expdata$Group, expdata["Donor"])
    dge <- cts %>% DGEList %>% calcNormFactors %>% estimateDisp(design=design, robust=TRUE)
    quantify.by.group.DGEList(dge, expdata$Group, expdata["Donor"])
})
s2.5.sampleFPKM.by.stype %<>% lapply(. %>% data.frame %>% {
    cbind(mcols(sexp2.5)[c("ENTREZ", "SYMBOL", "GENENAME")] %>% as("data.frame"), .)
})
s2.5.groupFPKM.by.stype %<>% lapply(. %>% data.frame %>% {
    cbind(mcols(sexp2.5)[c("ENTREZ", "SYMBOL", "GENENAME")] %>% as("data.frame"), .)
})

sample.fpkm.results <-
    list(H3K4me3=s1.sampleFPKM.by.stype$H3K4me3,
         H3K4me2=s1.sampleFPKM.by.stype$H3K4me2,
         H3K27me3=s2.5.sampleFPKM.by.stype$H3K27me3,
         input1KB=s1.sampleFPKM.by.stype$input,
         input2.5KB=s2.5.sampleFPKM.by.stype$input)
## This would be way too big probably
## write.xlsx.multisheet(sample.fpkm.results, "results/chipseq-sample-FPKM.xlsx")

group.fpkm.results <-
    list(H3K4me3=s1.groupFPKM.by.stype$H3K4me3,
         H3K4me2=s1.groupFPKM.by.stype$H3K4me2,
         H3K27me3=s2.5.groupFPKM.by.stype$H3K27me3,
         input1KB=s1.groupFPKM.by.stype$input,
         input2.5KB=s2.5.groupFPKM.by.stype$input)
write.xlsx.multisheet(group.fpkm.results, "results/chipseq-group-FPKM.xlsx")
