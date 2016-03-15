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

require(grid)

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

do.maplot <- function(tab) {
    tab %<>% arrange(PValue, logCPM.mean)
    ggplot(tab) + aes(x=logCPM.mean, y=logFC) +
        geom_hline(yintercept=0, alpha=0.5) +
        geom_point(aes(size=-log10(PValue), color=FDR <= 0.1)) +
        geom_density2d(show_guide=FALSE) +
        geom_smooth(color="red") +
        xlab("log2(CPM)") + ylab("log2(FC)") +
        scale_size(
            range=c(.1,1.5),
            trans="sqrt",
            name="-log10(PValue)") +
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

alltestnames <- c("Naive.AllT", "Memory.AllT", "Naive.D0vD1",
                  "Naive.D0vD5", "Naive.D0vD14", "Memory.D0vD1",
                  "Memory.D0vD5", "Memory.D0vD14", "NvM.AllT", "NvM.D0",
                  "NvM.D1", "NvM.D5", "NvM.D14", "Fac.AllT", "Fac.D1",
                  "Fac.D5", "Fac.D14", "MD0vND14")

{
    tsmsg("Loading 1kb promoter data")
    promoter.1kb.alltabs <- local({
        edger.analyses <- load.object.from.file("promoter-edger-results3-groups-1kb.RDa", "edger.analyses")
        tsmsg("Extracting promoter result tables")
        all.results.ql <- do.call(c, unname(plyr::llply(edger.analyses, `[[`, "results.ql")))
        rm(edger.analyses)
        gc()
        all.results.ql <- plyr::llply(all.results.ql, function(x) {
            x$Input.signif <- rownames(x) %in% rownames(all.results.ql$input.consistency)
            reorder.columns(x, start=c("ENTREZ", "SYMBOL", "GENENAME"))
        })
        all.results.ql[str_detect(names(all.results.ql), "H3K4me[23]")]
    })
    names(promoter.1kb.alltabs) %<>%
        str_replace("M0vN14", "MD0vND14") %>%
        str_replace("\\.vs\\.", "v") %>%
        str_replace("Fac\\.D0v", "Fac.")

    tsmsg("Loading 2.5kb promoter data")
    promoter.2.5kb.alltabs <- {
        load.object.from.file(
            "promoter-edger-results3-groups-2.5kb.RDa",
            "edger.analysis")$results.ql
    } %>% .[str_detect(names(.), "H3K27me3")]
    names(promoter.2.5kb.alltabs) %<>%
        str_replace("M0vN14", "MD0vND14") %>%
        str_replace("\\.vs\\.", "v") %>%
        str_replace("Fac\\.D0v", "Fac.")
    tsmsg("Loading RNA-seq tables")
    rnaseq.alltabs <- load.object.from.file("rnaseq-edger-results3.RDa", "edger.analysis")$results.ql


    tsmsg("Grouping tables by sample type")
    me2tabs <- llply(alltestnames, function(i) promoter.1kb.alltabs[[str_c("H3K4me2.", i)]])
    me3tabs <- llply(alltestnames, function(i) promoter.1kb.alltabs[[str_c("H3K4me3.", i)]])
    k27tabs <- llply(alltestnames, function(i) promoter.2.5kb.alltabs[[str_c("H3K27me3.", i)]])
    rnatabs <- rnaseq.alltabs[alltestnames]
    stopifnot(!any(sapply(list(me2tabs, me3tabs, k27tabs, rnatabs), is.null)))
    stopifnot(!any(sapply(c(me2tabs, me3tabs, k27tabs, rnatabs), is.null)))

    rm(promoter.1kb.alltabs, promoter.2.5kb.alltabs, rnaseq.alltabs); gc()
}

alltabs <-
    c(H3K4me2=me2tabs,
      H3K4me3=me3tabs,
      H3K27me3=k27tabs,
      RNASeq=rnatabs) %>%
    Filter(f=function(x) "logFC" %in% colnames(x))

{
    p <- lapply(names(alltabs), function(i) {
        do.maplot(alltabs[[i]]) +
            ggtitle(sprintf("MA plot for %s", i))
    })
    pdf("results/maplots.pdf", width=8, height=8)
    suppressMessages(print(p))
    dev.off()
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
           Sample=interaction(Sampletype, Group, Donor)) %>%
    droplevels

colData(sexp.1kb) %<>% fix.cd
colData(sexp.2.5kb) %<>% fix.cd
sexp.1kb %<>% set_colnames(colData(.)$Sample %>% as.character %>% make.unique) %>%
    .[,colData(.) %$% order(Sampletype, Celltype, Timepoint, Donor)]
sexp.2.5kb %<>% set_colnames(colData(.)$Sample %>% as.character %>% make.unique) %>%
    .[,colData(.) %$% order(Sampletype, Celltype, Timepoint, Donor)]

{
    individual.plots <- list()
    group.plots <- list()
    for (i in levels(colData(sexp.1kb)$Sampletype)) {
        sexp.sub <- sexp.1kb[,colData(sexp.1kb)$Sampletype == i]
        group <- colData(sexp.sub)$Group
        group.indices <- split(seq_len(ncol(sexp.sub)), group)
        tsmsg("Creating DGEList for ", i)
        dge <- DGEList(assay(sexp.sub)) %>%
            calcNormFactors %>%
            .[aveLogCPM(.) >= 1.5,] %>%
            estimateDisp(design=model.matrix(~Group + Donor, data=colData(sexp.sub)))
        logCPM.mean <- aveLogCPM(dge, prior.count=2)
        logCPM.samples <- cpm(dge, log=TRUE, prior.count=2)
        logCPM.group <-
            lapply(group.indices,
                   . %>% dge[,.] %>% aveLogCPM(prior.count=2)) %>%
            do.call(what=cbind)
        tsmsg("Creating plots for ", i)
        individual.plots %<>% c(lapply(colnames(logCPM.samples), function(s) {
            do.sample.maplot(logCPM.samples, s, logCPM.mean) +
                ggtitle(sprintf("MA plot for %s vs mean of other %s samples",
                                s, i))
        }))
        group.plots %<>% c(lapply(colnames(logCPM.group), function(s) {
            do.sample.maplot(logCPM.group, s, logCPM.mean) +
                ggtitle(sprintf("MA plot for %s vs mean of other %s group",
                                s, i))
        }))
    }
    for (i in levels(colData(sexp.2.5kb)$Sampletype)) {
        sexp.sub <- sexp.2.5kb[,colData(sexp.2.5kb)$Sampletype == i]
        group <- colData(sexp.sub)$Group
        group.indices <- split(seq_len(ncol(sexp.sub)), group)
        tsmsg("Creating DGEList for ", i)
        dge <- DGEList(assay(sexp.sub)) %>%
            calcNormFactors %>%
            .[aveLogCPM(.) >= 1.5,] %>%
            estimateDisp(design=model.matrix(~Group + Donor, data=colData(sexp.sub)))
        logCPM.mean <- aveLogCPM(dge, prior.count=2)
        logCPM.samples <- cpm(dge, log=TRUE, prior.count=2)
        logCPM.group <-
            lapply(group.indices,
                   . %>% dge[,.] %>% aveLogCPM(prior.count=2)) %>%
            do.call(what=cbind)
        tsmsg("Creating plots for ", i)
        individual.plots %<>% c(lapply(colnames(logCPM.samples), function(s) {
            do.sample.maplot(logCPM.samples, s, logCPM.mean) +
                ggtitle(sprintf("MA plot for %s vs mean of other %s samples",
                                s, i))
        }))
        group.plots %<>% c(lapply(colnames(logCPM.group), function(s) {
            do.sample.maplot(logCPM.group, s, logCPM.mean) +
                ggtitle(sprintf("MA plot for %s vs mean of other %s groups",
                                s, i))
        }))
    }
    tsmsg("Saving plots")
    pdf("results/sample-maplots.pdf", width=8, height=8)
    suppressMessages(print(individual.plots))
    dev.off()
    pdf("results/group-maplots.pdf", width=8, height=8)
    suppressMessages(print(group.plots))
    dev.off()
}

{
    system2("convert", c("-density", "600", "results/maplots.pdf", "results/maplots-raster.pdf"))
    system2("convert", c("-density", "600", "results/sample-maplots.pdf", "results/sample-maplots-raster.pdf"))
    system2("convert", c("-density", "600", "results/group-maplots.pdf", "results/group-maplots-raster.pdf"))
}
