#!/usr/bin/Rscript

source("rnaseq-common2.R", chdir=TRUE)

library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(annotate)
library(rtracklayer)
library(magrittr)
library(dplyr)
library(DESeq2)
require(grid)
library(scales)

flatten.CharacterList <- function(x, sep=",") {
    tryCatch(x <- as(x, "CharacterList"),
             error=function(...) NULL)
    if (is(x, "CharacterList"))
        setNames(rtracklayer:::pasteCollapse(x, collapse=sep),
                 names(x))
    else
        x
}

fix.cl.in.df <- function(df) {
    df[] <- llply(df, flatten.CharacterList, .parallel=TRUE)
}

load.object.from.file <- function(datafile, varname) {
    within(list(), load(datafile))[[varname]]
}


captureToRaster <- function(expr, ...) {
    dev.new(...)
    curdev <- dev.cur()
    on.exit(dev.off(curdev))
    expr
    grid.cap()
}

ggToRaster <- function(ggp, ...) {
    captureToRaster(print(ggp), ...)
}

## Minor breaks at every 1 logFC
maplot.minor.break.fun <- function(limits) {
    seq(from=ceiling(min(limits)),
        to=floor(max(limits)),
        by=1)
}

mutate_.DataFrame <- function(.data, ..., .dots) {
    .data %>% as.data.frame %>% mutate_(.data=., ..., .dots=.dots) %>% as("DataFrame")
}

## Now do per-sample MA plots
do.sample.maplot <- function(logCPM, samp, logCPM.mean=rowMeans(logCPM)) {
    stopifnot(ncol(logCPM) >= 2)
    if (is.logical(samp)) {
        samp %<>% which
    } else  if (is.character(samp)) {
        samp %<>% match(colnames(logCPM))
    }
    stopifnot(is.numeric(samp),
              length(samp) == 1,
              samp > 0,
              samp <= ncol(logCPM))

    tab <- data.frame(Samp=logCPM[,samp],
                      Other=rowMeans(logCPM[,-samp]),
                      logCPM.mean=logCPM.mean) %>%
        mutate(logFC=Samp - Other) %>%
        arrange(logCPM.mean)
    ggplot(tab) + aes(x=logCPM.mean, y=logFC) +
        geom_hline(yintercept=0, alpha=0.5) +
        geom_point(size=0.5) +
        geom_density2d(show_guide=FALSE) +
        geom_smooth(color="red") +
        xlab("log2(CPM)") + ylab("log2(FC)") +
        scale_x_continuous(
            expand=c(0.05, 0),
            ## labels=function(x) parse(text=sprintf("2^%s", x)),
            minor_breaks=maplot.minor.break.fun) +
        scale_y_continuous(
            expand=c(0.05, 0),
            ## labels=function(x) parse(text=sprintf("2^%s", x)),
            minor_breaks=maplot.minor.break.fun,
            limits=c(-1, 1) * max(abs(tab$logFC))) +
        theme_bw() +
        theme(
            legend.position=c(1,1),
            legend.justification=c(1.1, 1.1))
}

sexp.1kb <- readRDS("promoter-group-counts-1kb.RDS") %>%
    updateObject %>%
    .[,colData(.)$Sampletype %in% c("H3K4me2", "H3K4me3")]
sexp.2.5kb <- readRDS("promoter-group-counts-2.5kb.RDS") %>%
    updateObject %>%
    .[,colData(.)$Sampletype == "H3K27me3"]

fix.cd <- . %>%
    mutate(
        Celltype=factor(Celltype, levels=c("Naive", "Memory")),
        Timepoint=factor(Timepoint, levels=c("D0", "D1", "D5", "D14")),
        Group=interaction(Celltype, Timepoint),
        Sample=interaction(Sampletype, Group, Donor)) %>% {
            levels(.$Timepoint) %<>% str_replace("^D", "Day ")
            levels(.$Group) %<>% str_replace("\\.D", " Day ")
            levels(.$Donor) %<>% str_replace("^([0-9]+)", "Donor \\1")
            levels(.$Celltype) %<>% str_replace("(Naive|Memory)$", "\\1 CD4 T-cells")
            .
        } %>%
    droplevels

colData(sexp.1kb) %<>% fix.cd
colData(sexp.2.5kb) %<>% fix.cd
sexp.1kb %<>% set_colnames(colData(.)$Sample %>% as.character %>% make.unique) %>%
    .[,colData(.) %$% order(Sampletype, Celltype, Timepoint, Donor)]
sexp.2.5kb %<>% set_colnames(colData(.)$Sample %>% as.character %>% make.unique) %>%
    .[,colData(.) %$% order(Sampletype, Celltype, Timepoint, Donor)]

sexps <- list(H3K4me2=sexp.1kb %>% .[,colData(.)$Sampletype == "H3K4me2"],
              H3K4me3=sexp.1kb %>% .[,colData(.)$Sampletype == "H3K4me3"],
              H3K27me3=sexp.2.5kb %>% .[,colData(.)$Sampletype == "H3K27me3"])
dgelists <- lapply(sexps, function(sexp) {
    dge <- DGEList(assay(sexp, "counts")) %>%
        calcNormFactors %>%
        .[aveLogCPM(.) >= 1.5,] %>%
        estimateDisp(design=model.matrix(~Group + Donor, data=colData(sexp)))
    dge$samples %<>% cbind(colData(sexp)) %>% as.data.frame
    dge
})

sample.madata <- lapply(dgelists, function(dge) {
    logCPM.mean <- aveLogCPM(dge, prior.count=2)
    logCPM.samples <- cpm(dge, log=TRUE, prior.count=2)
    lapply(seq_len(ncol(logCPM.samples)), function(i) {
        samp <- colnames(logCPM.samples)[i]
        tsmsg("Generating data table for ", samp)
        data.frame(SampleID=samp,
                   SampleLogCPM=logCPM.samples[,i],
                   OtherLogCPM=aveLogCPM(dge[,-i], prior.count=2),
                   MeanLogCPM=logCPM.mean) %>%
            mutate(logFC=SampleLogCPM - OtherLogCPM) %>%
            cbind(dge$samples[i,] %>% set_rownames(NULL)) %>%
            arrange(MeanLogCPM)
    }) %>% do.call(what=rbind)
}) %>% do.call(what=rbind)

p32 <- sample.madata %>% group_by(Sampletype) %>% do(Plot={
    df <- .
    ggplot(df) + aes(x=MeanLogCPM, y=logFC) +
        geom_hline(yintercept=0, alpha=0.5) +
        geom_point(size=0.5) +
        geom_density2d(show_guide=FALSE) +
        geom_smooth(color="red") +
        xlab("log2(CPM)") + ylab("log2(FC)") +
        coord_equal(
            xlim=c(min(df$MeanLogCPM), 8),
            ylim=c(-1, 1) * quantile(abs(df$logFC), 0.95)) +
        scale_x_continuous(
            minor_breaks=maplot.minor.break.fun) +
        scale_y_continuous(
            minor_breaks=maplot.minor.break.fun) +
        theme_bw() +
        facet_grid(Donor ~ Group) +
        ggtitle(sprintf("MA plots for %s samples", df$Sampletype[1]))
}, SmoothPlot={
    df <- .
    ggplot(df) + aes(x=MeanLogCPM, y=logFC) +
        ## Smooth scatter
        stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
        scale_fill_gradient(low="white", high=muted("blue")) +
        ## Contour plot
        stat_density2d(colour="black", size=.25, show_guide=FALSE) +
        geom_hline(yintercept=0, alpha=0.5) +
        geom_smooth(color="red") +
        xlab("log2(CPM)") + ylab("log2(FC)") +
        coord_equal(
            xlim=c(min(df$MeanLogCPM), 8),
            ylim=c(-1, 1) * quantile(abs(df$logFC), 0.95)) +
        scale_x_continuous(
            minor_breaks=maplot.minor.break.fun) +
        scale_y_continuous(
            minor_breaks=maplot.minor.break.fun) +
        theme_bw() +
        facet_grid(Donor ~ Group) +
        ggtitle(sprintf("MA plots for %s samples", df$Sampletype[1]))
})
{
    pdf("results/sample-maplots-grid32.pdf", width=20, height=8)
    print(p32$Plot)
    dev.off()
}
{
    pdf("results/sample-maplots-grid32-smooth.pdf", width=21, height=8)
    print(p32$SmoothPlot)
    dev.off()
}

p16 <- sample.madata %>% group_by(Sampletype, Celltype) %>% do(Plot={
    df <- .
    ggplot(df) + aes(x=MeanLogCPM, y=logFC) +
        geom_hline(yintercept=0, alpha=0.5) +
        geom_point(size=0.5) +
        geom_density2d(show_guide=FALSE) +
        geom_smooth(color="red") +
        xlab("log2(CPM)") + ylab("log2(FC)") +
        coord_equal(
            xlim=c(min(df$MeanLogCPM), 8),
            ylim=c(-1, 1) * quantile(abs(df$logFC), 0.95)) +
        scale_x_continuous(
            minor_breaks=maplot.minor.break.fun) +
        scale_y_continuous(
            minor_breaks=maplot.minor.break.fun) +
        theme_bw() +
        facet_grid(Donor ~ Timepoint) +
        ggtitle(sprintf("MA plots for %s samples in %s", df$Sampletype[1], df$Celltype[1]))
}, SmoothPlot={
    df <- .
    ggplot(df) + aes(x=MeanLogCPM, y=logFC) +
        ## Smooth scatter
        stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
        scale_fill_gradient(low="white", high=muted("blue")) +
        ## Contour plot
        stat_density2d(colour="black", size=.25, show_guide=FALSE) +
        geom_hline(yintercept=0, alpha=0.5) +
        geom_smooth(color="red") +
        xlab("log2(CPM)") + ylab("log2(FC)") +
        coord_equal(
            xlim=c(min(df$MeanLogCPM), 8),
            ylim=c(-1, 1) * quantile(abs(df$logFC), 0.95)) +
        scale_x_continuous(
            minor_breaks=maplot.minor.break.fun) +
        scale_y_continuous(
            minor_breaks=maplot.minor.break.fun) +
        theme_bw() +
        facet_grid(Donor ~ Timepoint) +
        ggtitle(sprintf("MA plots for %s samples in %s", df$Sampletype[1], df$Celltype[1]))
})
{
    pdf("results/sample-maplots-grid16.pdf", width=10, height=8)
    print(p16$Plot)
    dev.off()
}
{
    pdf("results/sample-maplots-grid16-smooth.pdf", width=11, height=8)
    print(p16$SmoothPlot)
    dev.off()
}

files <- c("results/sample-maplots-grid32.pdf",
           "results/sample-maplots-grid16.pdf")
for (infile in files) {
    outfile <- infile %>% str_replace("\\.pdf$", "-raster.pdf")
    tsmsg("Converting ", infile, ", to raster")
    system2("convert", c("-density", "600", infile, outfile))
}
