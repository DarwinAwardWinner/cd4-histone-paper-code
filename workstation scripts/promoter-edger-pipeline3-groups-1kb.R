#!/usr/bin/Rscript

library(Biobase)
library(RColorBrewer)
library(Matrix)
library(rtracklayer)
library(GenomicRanges)
library(magrittr)
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

sexp <- readRDS("promoter-group-counts-1kb.RDS")
## Load CPGi data and compute overlaps
cpgi <- read.table("CpGi.Ext.bed", sep="\t")
cpgi <- GRanges(seqnames=cpgi[[1]], ranges=IRanges(start=cpgi[[2]]+1, end=cpgi[[3]]), strand="*")
mcols(sexp)$CpGi <- overlapsAny(sexp, cpgi)

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

## MDS Plots (techincal replicates)
{
    tsmsg("Generating MDS plots for techincal replicates")
    pdf("results/1kb/MDSplots_techrep.pdf")
    tryCatch({
        dge <- DGEList(counts)
        pal <- brewer.pal(9, "Set1")
        always.label <- c("Celltype", "Sampletype", "Timepoint")
        mds <- NULL
        for (i in c("Run", "Donor", "Celltype", "Sampletype", "Timepoint")) {
            labels <- with(expdata, Reduce(`:`, mget(union(i, always.label))))
            tsmsg("Plotting MDS with colors = ", i)
            mds <- plotMDS((if (is.null(mds)) dge else mds), labels=labels, main=sprintf("MDS Plot (Colors = %s)", i),
                           col=pal[as.numeric(as.factor(expdata[[i]]))],
                           cex=0.5)
        }
    }, finally=dev.off())
    tsmsg("Finished generating MDS plots")
}

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
    x <- droplevels(Sampletype:Celltype:Timepoint)
    ## Combine all input samples into one group
    ## x[Sampletype == "input"] <- "input"
    x
})

## MDS Plots (biological replicates)
{
    tsmsg("Generating MDS plots")
    pal <- brewer.pal(9, "Set1")
    always.label <- c("Celltype", "Sampletype", "Timepoint")
    dge <- DGEList(counts)
    pdf("results/1kb/MDSplots_biolrep.pdf")
    tryCatch({
        mds <- NULL
        for (i in c("Donor", "Celltype", "Sampletype", "Timepoint")) {
            labels <- with(expdata, Reduce(`:`, mget(union(i, always.label))))
            tsmsg("Plotting MDS with colors = ", i)
            mds <- plotMDS((if (is.null(mds)) dge else mds), labels=labels, main=sprintf("MDS Plot (Colors = %s)", i),
                           col=pal[as.numeric(as.factor(expdata[[i]]))],
                           cex=0.5)
        }
    }, finally=dev.off())
    outliers <- c("4659.Naive.H3K4me3.D1",
                  "5291.Memory.H3K27me3.D14",
                  "4659.Naive.input.D1",
                  "5291.Memory.input.D1",
                  "5291.Naive.input.D1",
                  "4659.Memory.input.D1",
                  "5291.Naive.input.D0")
    llply(levels(expdata$Sampletype), function(stype) {
        sel <- expdata$Sampletype == stype & ! rownames(expdata) %in% outliers
        always.label <- c("Celltype", "Timepoint")
        dge <- DGEList(counts[,sel])
        mds <- NULL
        pdf(sprintf("results/1kb/MDSplots_%s_only.pdf", stype))
        tryCatch({
            for (i in c("Donor", "Celltype", "Timepoint")) {
                labels <- with(expdata[sel,], Reduce(`:`, mget(union(i, always.label))))
                tsmsg(sprintf("Plotting MDS for %s with colors = %s", stype, i))
                mds <- plotMDS((if (is.null(mds)) dge else mds), labels=labels, main=sprintf("MDS Plot (Colors = %s)", i),
                               col=pal[as.numeric(as.factor(expdata[[i]][sel]))],
                               cex=0.5)
            }
        }, finally=dev.off())
    }, .parallel=TRUE)
    tsmsg("Finished generating MDS plots")
}

## Verification of consistent dispersion estimates across "Donor",
## "Celltype", "Sampletype", "Timepoint"
{
    tsmsg("Generating factor dispersion plots")
    ## rownames(design) <- rownames(expdata)
    group.by.factor <- function(f) split(seq_along(f), f)

    groups <-
        with(expdata,
             c(list(All=seq(nrow(expdata))),
               Donor=group.by.factor(Donor),
               Celltype=group.by.factor(Celltype),
               Sampletype=group.by.factor(Sampletype),
               Timepoint=group.by.factor(Timepoint)))
    ## group.designs <- llply(groups, function (g) {
    ##     subdesign(design, g)
    ## })

    group.dgelists <- llply(names(groups), function(i) {
        ## gdesign <- group.designs[[i]]
        tsmsg("Processing group ", i)
        dge <- DGEList(counts[,groups[[i]]])
        tsmsg("Calculating norm factors for group ", i)
        dge <- calcNormFactors(dge)
        tsmsg("Estimating dispersions for group ", i)
        dge <- estimateDisp(dge, prior.df=0, robust=TRUE)
        tsmsg("Common disp for group ", i, ": ", dge$common.dispersion)
        dge
    }, .parallel=TRUE)
    names(group.dgelists) <- names(groups)

    group.bcv <- sqrt(sapply(group.dgelists, `[[`, "common.dispersion"))
    group.meancount <- sapply(group.dgelists, function(x) mean(x$samples$lib.size))
    group.bcv.trends <- llply(names(group.dgelists), function(gname) {
        dge <- group.dgelists[[gname]]
        A <- dge$AveLogCPM
        tdisp <- dge$trended.dispersion
        unique(data.frame(logCPM=A,
                          BCV=sqrt(tdisp),
                          Group=factor(gname, levels=names(group.dgelists))))
    }, .parallel=TRUE)
    all.bcv.trends <- do.call(rbind, group.bcv.trends)
    factor.bcv.trends <- list()
    for (fac in c("Donor", "Celltype", "Sampletype", "Timepoint")) {
        sel <- str_detect(names(group.bcv.trends), perl(fac))
        factor.bcv.trends[[fac]] <- do.call(rbind, group.bcv.trends[sel])
    }

    pdf("results/1kb/Factor BCV trends.pdf")
    print(ggplot(all.bcv.trends) + aes(x=logCPM, y=BCV, group=Group, colour=Group) + geom_line() +
          ggtitle("Factor BCV Trends"))
    for (fac in names(factor.bcv.trends)) {
        print(ggplot(factor.bcv.trends[[fac]]) + aes(x=logCPM, y=BCV, group=Group, colour=Group) + geom_line() +
              ggtitle(sprintf("%s BCV Trends", fac)))
    }
    dev.off()
    write.xlsx(data.frame(BCV=group.bcv, Zscore=(group.bcv-group.bcv["All"])/sd(group.bcv[-1])),
               "results/1kb/Factor BCV.xlsx")
    tsmsg("Finished generating factor dispersion plots")
}

## Verification of consistent dispersion estimates between relevant groups
{
    tsmsg("Generating group dispersion plots")
    dge <- DGEList(counts)

    groups <- group.by.factor(expdata$Group)

    group.dgelists <- llply(names(groups), function(i) {
        tsmsg("Processing group ", i)
        dge <- dge[,groups[[i]]]
        tsmsg("Calculating norm factors for group ", i)
        dge <- calcNormFactors(dge)
        tsmsg("Estimating dispersions for group ", i)
        repeat {
            proc <- mcparallel(estimateDisp(dge))
            res <- mccollect(proc)[[1]]
            if (is.null(res$common.dispersion))
                message("Dispersion estimation crashed. Retrying.")
            else {
                dge <- res
                break
            }
        }
        ## dge <- estimateDisp(dge)
        tsmsg("Common disp for group ", i, ": ", dge$common.dispersion)
        dge
    }, .parallel=TRUE)
    names(group.dgelists) <- names(groups)

    group.bcv <- sqrt(sapply(group.dgelists, `[[`, "common.dispersion"))
    group.bcv.trends <- lapply(names(group.dgelists), function(gname) {
        dge <- group.dgelists[[gname]]
        A <- dge$AveLogCPM
        tdisp <- dge$trended.dispersion
        unique(data.frame(logCPM=A,
                          BCV=sqrt(tdisp),
                          Group=factor(gname, levels=names(group.dgelists))))
    })
    all.bcv.trends <- do.call(rbind, group.bcv.trends)
    all.facs <- c("Sampletype", "Timepoint", "Celltype")
    all.bcv.trends <- data.frame(all.bcv.trends, expdata[match(all.bcv.trends$Group, expdata$Group),
                                                         all.facs])
    pdf("results/1kb/Group BCV trends.pdf")
    tsmsg("Printing all-group BCV plot")
    print(ggplot(all.bcv.trends) +
          aes(x=logCPM, y=BCV, group=Group, colour=Sampletype:Celltype, linetype=Timepoint) +
          geom_line() +
          ggtitle("All Group BCV Trends"))
    for (fac in all.facs) {
        for (val in levels(all.bcv.trends[[fac]])) {
            tsmsgf("Printing BCV plot for %s = %s", fac, val)
            data <- droplevels(all.bcv.trends[all.bcv.trends[[fac]] == val,])
            otherfacs <- setdiff(all.facs, fac)
            data$color <- data[[otherfacs[1]]]
            data$linetype <- data[[otherfacs[2]]]
            print(ggplot(data) +
                  aes(x=logCPM, y=BCV, group=Group, colour=color, linetype=linetype) + geom_line() +
                  ggtitle(sprintf("Group BCV Trends (%s = %s)", fac, val)))
        }
    }
    dev.off()
    group.bcv.table <- data.frame(BCV=group.bcv, expdata[match(names(group.bcv), expdata$Group), all.facs])
    factor.bcv.summaries <-
        with(group.bcv.table,
             llply(all.facs, function(fac) {
                 data.frame(sapply(split(BCV, get(fac)), summary))
             }))
    names(factor.bcv.summaries) <- all.facs
    fit.with.input <- lm(BCV ~ Sampletype + Timepoint + Celltype,
                         data=group.bcv.table)
    test.with.input <- anova(fit.with.input)
    fit <- lm(BCV ~ Sampletype + Timepoint + Celltype,
              data=droplevels(group.bcv.table[group.bcv.table$Sampletype != "input",]))
    test <- anova(fit)
    write.xlsx.multisheet(c(list(`All Groups`=group.bcv.table),
                            factor.bcv.summaries,
                            list(ANOVA=as.data.frame(test),
                                 `ANOVA with input`=as.data.frame(test.with.input))),
                          "results/1kb/Group BCV.xlsx")
    tsmsg("Finished generating group dispersion plots")
}

## Differential binding tests
{
    design <- model.matrix(terms(~0 + Timepoint:Celltype:Sampletype + Sampletype:Donor, keep.order=TRUE), droplevels(expdata))
    colnames(design)[1:nlevels(expdata$Group)] <- levels(droplevels(expdata$Group))
    colnames(design) <- make.names(colnames(design))

    ## design2 <- make.group.design(expdata$Group, block=expdata["Donor"])
    input.groups <- grep("^input\\.", colnames(design), value=TRUE)
    all.timepoints <- levels(expdata$Timepoint)
    nonzero.timepoints <- setdiff(all.timepoints, "D0")
    ip.names <- setdiff(levels(expdata$Sampletype), "input")

    alltests <- list()
    alltests$input.consistency <- all.contrasts(input.groups)
    for (ip in c("H3K4me2", "H3K4me3")) {
        ## Timepoint tests
        for (ctype in c("Naive", "Memory")) {
            coefs <- sprintf("%s.%s.%s", ip, ctype, all.timepoints)
            prefix <- sprintf("%s.%s.", ip, ctype)
            ## ANOVA tests for all timepoints
            alltests[[sprintf("%s.%s.AllT", ip, ctype)]] <- {
                all.contrasts(coefs, prefix)
            }
            ## Pairwise tests for every pair of timepoints
            alltests %<>% c(all.pairwise.contrasts(coefs, prefix))
        }
        1; {
            ## NvM tests for all/each timepoint
            nvm.ct <- sprintf("%s.Memory.%s - %s.Naive.%s", ip, all.timepoints, ip, all.timepoints) %>%
                setNames(sprintf("%s.NvM.%s", ip, all.timepoints))
            ## ANOVA
            alltests[[sprintf("%s.NvM.AllT", ip)]] <- nvm.ct
            ## Individual
            alltests %<>% c(nvm.ct)
        }
        ## Factorial tests
        1; {
            ## ANOVA of all factorial contrasts
            fact.ct <- all.contrasts(nvm.ct)
            names(fact.ct) <- sprintf("%s.Fac.%s", ip, names(all.contrasts(all.timepoints, NULL)))
            alltests[[sprintf("%s.Fac.AllT", ip)]] <- fact.ct

            ## Factorial tests for each pair of timepoints
            fact.pw.ct <- all.pairwise.contrasts(nvm.ct)
            names(fact.pw.ct) <- sprintf("%s.Fac.%s", ip, names(all.pairwise.contrasts(all.timepoints, NULL)))
            alltests %<>% c(fact.pw.ct)
        }
        ## Naive D14 vs Mem D0
        alltests[[sprintf("%s.MD0vND14", ip)]] <- sprintf("%s.Memory.D0 - %s.Naive.D14", ip, ip)
        ## Overall test for diffexp
        alltests[[sprintf("%s.AllExp", ip)]] <-
            expand.grid(ctype=c("Naive", "Memory"), tp=all.timepoints) %>%
            mlply(function(ctype, tp) {
                sprintf("%s.%s.%s", ip, ctype, tp)
            }) %>% unlist %>% all.contrasts
    }

    donor.tests <-
        setNames(lapply(levels(expdata$Sampletype), function(stype) {
            colnames(design)[str_detect(colnames(design), fixed(str_c(stype, ".Donor")))]
        }), str_c(levels(expdata$Sampletype), ".Donor"))

    alltests %<>% c(donor.tests)

    included.timepoints <- llply(alltests, function(test) {
        if (any(str_detect(test, "Donor"))) {
            ## Donor tests should be conducted on all timepoints
            sort(levels(expdata$Timepoint))
        } else {
            ## Other tests should be conducted only on the timepoints
            ## involved, since dispersion varies widely across
            ## timepoints.
            sort(unique(unlist(str_extract_all(test, perl("D\\d+")))))
        }
    })

    testnames.by.tpgroup <- split(names(included.timepoints), sapply(included.timepoints, str_c, collapse="."))
    timepoints.by.tpgroup <- setNames(included.timepoints[sapply(testnames.by.tpgroup, `[`, 1)],
                                      names(testnames.by.tpgroup))
    alltests.by.tpgroup <- llply(testnames.by.tpgroup, function(i) alltests[i])
    stopifnot(all(sort(unique(unlist(included.timepoints))) %in% levels(expdata$Timepoint)))
    stopifnot(min(str_length(names(testnames.by.tpgroup))) > 0)

    ## Want at least 5 counts in at least one group (4 samples)
    groupsize <- 4
    feature.subset <- function(counts, ...) {
        rowQ(counts, ncol(counts) - groupsize + 1) >= 5
    }

    edger.analyses <-
        llply(names(testnames.by.tpgroup), function(tpgroup.name) {
            incl.tp <- timepoints.by.tpgroup[[tpgroup.name]]
            incl.tests <- alltests.by.tpgroup[[tpgroup.name]]
            tsmsgf("Doing tests for Timepoint group %s", deparse(incl.tp))
            selected.samples <- expdata$Timepoint %in% incl.tp
            tsmsgf("Selected %s out of %s samples for this group",
                   sum(selected.samples), ncol(counts))
            subdes <- subdesign(design, selected.samples, intcol=NULL)
            stopifnot(min(dim(subdes)) > 0)
            x <- do.edger.analysis(counts[,selected.samples], subdes, incl.tests,
                                   feature.annot=featuredata,
                                   feature.subset=feature.subset,
                                   group=expdata$Group[selected.samples],
                                   .parallel=FALSE)
            tsmsgf("Finished tests for Timepoint group %s", deparse(incl.tp))
            x[c("dge", "results.ql", "results.lr", "topresults.ql", "topresults.lr")]
        }, .parallel=FALSE)
    names(edger.analyses) <- names(testnames.by.tpgroup)

    tsmsg("Extracting result tables")
    all.topresults.ql <- do.call(c, unname(llply(edger.analyses, `[[`, "topresults.ql")))[names(alltests)]
    all.topresults.ql <- llply(all.topresults.ql, function(x) {
        x$Input.signif <- rownames(x) %in% rownames(all.topresults.ql$input.consistency)
        reorder.columns(x, start=c("ENTREZ", "SYMBOL", "GENENAME"))
    })

    all.topresults.lr <- do.call(c, unname(llply(edger.analyses, `[[`, "topresults.lr")))[names(alltests)]
    all.topresults.lr <- llply(all.topresults.lr, function(x) {
        ## More stringent cutoff
        x <- x[x$FDR <= 0.05,]
        x$Input.signif <- rownames(x) %in% rownames(all.topresults.lr$input.consistency)
        reorder.columns(x, start=c("ENTREZ", "SYMBOL", "GENENAME"))
    })

    diffcounts.ql <- data.frame(cbind(diffcount=sapply(all.topresults.ql, nrow)))
    diffcounts.lr <- data.frame(cbind(diffcount=sapply(all.topresults.lr, nrow)))

    tsmsg("Saving R data")
    save.image("promoter-edger-results3-groups-1kb.RDa")

    tsmsg("Saving result tables")
    write.xlsx.multisheet(addsummary(all.topresults.ql),
                          "results/1kb/promoter-edger-topgroups3-ql.xlsx",
                          row.names=TRUE)
    ## TOO MANY GENES! Crashes Java!
    ## try({
    ##     write.xlsx.multisheet(addsummary(all.topresults.lr),
    ##                           "results/1kb/promoter-edger-topgroups3-lr.xlsx",
    ##                           row.names=TRUE)
    ## })



    tsmsg("Plotting p-value distributions")
    binbreaks <- seq(0, 1, length.out=41)
    smallbinbreaks <- seq(0, .1, length.out=5)
    binbreaks <- seq(0, 1, by=1/80)
    bin.df <- data_frame(
        leftedge = binbreaks[-length(binbreaks)],
        rightedge = binbreaks[-1],
        width = rightedge - leftedge,
        midpoint = (rightedge + leftedge) / 2)

    all.results.lr <- edger.analyses %>% lapply(`[[`, "results.lr") %>%
        unname %>% do.call(what=c) %>% .[names(alltests)]
    all.results.ql <- edger.analyses %>% lapply(`[[`, "results.ql") %>%
        unname %>% do.call(what=c) %>% .[names(alltests)]
    bin.tables.lr <- all.results.lr %>%
        lapply(. %>%
               filter(logCPM.mean >= 1) %$%
               PValue %>%
               cut(breaks=binbreaks) %>%
               table %>% as.data.frame %>%
               setNames(c("Bin", "PCount")) %>%
               cbind(bin.df) %>%
               mutate(PDensity=PCount / sum(PCount) / width))
    bin.tables.ql <- all.results.ql %>%
        lapply(. %>%
               filter(logCPM.mean >= 1) %$%
               PValue %>%
               cut(breaks=binbreaks) %>%
               table %>% as.data.frame %>%
               setNames(c("Bin", "PCount")) %>%
               cbind(bin.df) %>%
               mutate(PDensity=PCount / sum(PCount) / width))

    pdf("results/1kb/Pvalue distributions LR.pdf")
    for (i in names(all.results.lr)) {
        print(ggplot(bin.tables.lr[[i]]) +
              aes(x=midpoint, y=PDensity, width=width) +
              geom_bar(stat="identity") +
              ## geom_histogram(aes(y = ..density..), breaks=binbreaks) +
              geom_hline(yintercept=1, alpha=0.5, linetype="dotted") +
              xlim(0,1) + ggtitle(sprintf("P-value distribution for %s test", i)) +
              xlab("p-value") + ylab("Relative frequency"))
    }
    dev.off()
    pdf("results/1kb/Pvalue distributions QL.pdf")
    for (i in names(all.results.ql)) {
        print(ggplot(bin.tables.ql[[i]]) +
              aes(x=midpoint, y=PDensity, width=width) +
              geom_bar(stat="identity") +
              ## geom_histogram(aes(y = ..density..), breaks=binbreaks) +
              geom_hline(yintercept=1, alpha=0.5, linetype="dotted") +
              xlim(0,1) + ggtitle(sprintf("P-value distribution for %s test", i)) +
              xlab("p-value") + ylab("Relative frequency"))
    }
    dev.off()
    pdf("results/1kb/PValue dist ExpVsDonor QL.pdf")
    for (ip in c("H3K4me2", "H3K4me3")) {
        i <- ip %>% str_c(c(".AllExp", ".Donor")) %>% setNames(c("Experimental Effects", "Donor Effects"))
        tb <- ldply(i, .id = "Test", . %>% {
            bin.tables.ql[[.]]
        }) %>%
            rbind(data.frame(
                Test="(Theoretical Uniform)",
                Bin=levels(.$Bin),
                PCount=1,
                bin.df,
                PDensity=1))
        print(ggplot(tb) +
              aes(x=midpoint, y=PDensity,
                  group=Test, linetype=Test) +
              geom_line(size=0.5) +
              scale_linetype_manual(values=c("solid", "dashed", "dotted")) +
              xlim(0,1) +
              ylim(0, max(tb$PDensity)) +
              ggtitle(sprintf("P-value distributions for %s", ip)) +
              xlab("p-value") + ylab("Relative frequency") +
              theme(legend.justification=c(1,1),
                    legend.position=c(0.95,0.95),
                    legend.background=element_rect(fill="gray90")))
    }
    dev.off()
}
