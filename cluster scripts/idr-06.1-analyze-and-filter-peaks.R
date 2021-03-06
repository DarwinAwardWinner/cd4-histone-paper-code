#!/bin/bash
# -*- mode:R -*-
#PBS
#PBS -j oe -o idr-06.1-analyze-and-filter-peaks.log
#SBATCH -c8
module load R/3.1.0;
cd $PBS_O_WORKDIR

Rscript --version
Rscript - <<'EOF'

library(annotate)
library(xlsx)
library(plyr)
options(cores=8)
options(mc.cores=8)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(Rsubread)
## Not sure which of these is the real option. Some code seems to
## use different ones.
options(cores=parallel:::detectCores())
options(mc.cores=parallel:::detectCores())
registerDoParallel()
library(GenomicRanges)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

## If llply is called on an unnamed character vector, use the
## character vector itself for the names.
llply <- function(.data, ...) {
    if (is.character(.data) && is.null(names(.data))) {
        names(.data) <- .data
    }
    plyr::llply(.data, ...)
}

fixdfchar <- function(df) {
    for (i in seq_along(df)) {
        if (is.factor(df[[i]]) && !any(duplicated(df[[i]]))) {
            df[[i]] <- as.character(df[[i]])
        }
    }
    df
}
splitargs <- function(...) {
    all.args <- unlist(list(...))
    unname(unlist(str_split(all.args, "\\s+")))
}

makecmd <- function(...) {
    args <- as.character(unlist(list(...)))
    cmd <- str_c(shQuote(args), collapse=" ")
}

whoami <- function() {
    system2("whoami", stdout=TRUE)
}

job.name.running <- function(name) {
    unlist(
        mclapply(name,
                 function(n)
                 length(system2("qselect", c("-N", n, "-s", "EHQRTW"), stdout=TRUE)) > 0,
                 mc.cores=50))
}

submit.job <- function(cmds, qsub.opts="", wd=getwd(),
                       name=str_replace_all(cmd, "[^A-Za-z0-9_-=+:,.!@#$%^&*()[]{}<>?|]", "_"),
                       logfile,
                       modules=c(), dry.run=FALSE,
                       allow.duplicate.name=FALSE) {
    if (!allow.duplicate.name && !is.null(name) && job.name.running(name)) {
        warning("Not starting duplicate job ", name)
        return()
    }
    if (!missing(logfile)) {
        qsub.opts <- c(qsub.opts, str_c("-j oe -o ", logfile))
    }
    temp <- tempfile()
    text <- c("#!/bin/sh",
              sprintf("#PBS -w %s -N %s", wd, name),
              str_c("#PBS ", str_c(qsub.opts, collapse=" ")),
              str_c("cd $PBS_O_WORKDIR"),
              sprintf("module load %s", modules),
              cmds)
    if (dry.run) {
        writeLines(text)
    }
    else {
        writeLines(text, con=temp)
        system2("qsub", args=temp)
        file.remove(temp)
    }
}

string.ends.with <- function(str, suffix) {
    str_sub(str, -str_length(suffix)) == suffix
}

file.name.as.directory <- function(f) {
    d <- file.path(f)
    if (! string.ends.with(d, .Platform$file.sep)) {
        d <- str_c(d, .Platform$file.sep)
    }
    d
}

data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE,
                               startOffset=0, endOffset=0) {
    stopifnot(class(df) == "data.frame")
    stopifnot(all(c("start", "end") %in% names(df)))
    stopifnot(any(c("chr", "seqnames") %in% names(df)))
    if("seqnames" %in% names(df))
        names(df)[names(df) == "seqnames"] <- "chr"
    if(!ignoreStrand && "strand" %in% names(df)) {
        if(is.numeric(df$strand)) {
            strand <- ifelse(df$strand == 1, "+", "*")
            strand[df$strand == -1] <- "-"
            df$strand <- strand
        }
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start + startOffset, end = df$end + endOffset),
                      strand = df$strand)
    } else {
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start + startOffset, end = df$end + endOffset))
    }
    if(keepColumns) {
        dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
                     "DataFrame")
        elementMetadata(gr) <- dt
    }
    names(gr) <- rownames(df)
    gr
}

read.narrowPeak <- function(file, ...) {
    peaks.df <- read.table(file, sep="\t", row.names=NULL, ...)
    names(peaks.df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
    peaks.df$strand <- "*"
    peaks.df$name <- as.character(peaks.df$name)
    row.names(peaks.df) <- str_replace(as.character(peaks.df$name), "^.*_peak", "peak")
    data.frame2GRanges(peaks.df, keepColumns=TRUE, startOffset=1, endOffset=0)
}

write.narrowPeak <- function(x, file, ...) {
    x <- as(x, "data.frame")
    if("seqnames" %in% names(x))
        names(x)[names(x) == "seqnames"] <- "chr"
    x <- x[c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")]
    write.table(x, file, sep="\t", row.names=FALSE, col.names=FALSE, ...)
}

wd <- "/gpfs/group/salomon/userdata/ryan/Projects/sarah-cd4"
setwd(wd)

qsub.opts <- "-l walltime=48:00:00"

peakcall.dir <- file.path(wd, "peak_calls")

idr.dir <- file.path(wd, "idr_analysis", "old")
dir.create(idr.dir, recursive=TRUE, showWarnings=FALSE)

bca.script <- "/gpfs/home/rcthomps/src/idrCode/batch-consistency-analysis.r"

tsmsg("Loading sample data")

expdata <- read.xlsx("sampledata.xlsx", 1)
expdata <- fixdfchar(expdata)
rownames(expdata) <- expdata$Sample
expdata$Lane <- factor(expdata$Lane)
expdata$Timepoint <- factor(expdata$Timepoint, levels=str_c("D", c(0, 1, 5, 14)))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
expdata$Sampletype <- factor(expdata$Sampletype, levels=c("input", "H3K4me2", "H3K4me3", "H3K27me3"))
expdata$Sample <-
    with(expdata, make.names(Sampletype:Celltype:Timepoint:Donor))
expdata$Group <- with(expdata, Sampletype:Celltype:Timepoint)

inputsamples <- expdata$Sampletype == "input"

expdata.noninput <- droplevels(expdata[!inputsamples,])
expdata.noninput$Peakfile <- file.path(peakcall.dir, sprintf("%s_peaks_sorted.encodePeak", expdata.noninput$Sample))

stopifnot(all(file.exists(expdata.noninput$Peakfile)))

idr.jobs <- ddply(expdata.noninput, .(Group), function(df) {
    stopifnot(nrow(df) >= 3)
    ## Get every pair of rows
    rowpairs <- combn(row.names(df), 2)
    Group <- df$Group[1]
    adply(rowpairs, 2, function(rn) {
        Donor1 <- df[rn[1],]$Donor
        Donor2 <- df[rn[2],]$Donor
        Sample1 <- df[rn[1],]$Sample
        Sample2 <- df[rn[2],]$Sample
        Peakfile1 <- df[rn[1],]$Peakfile
        Peakfile2 <- df[rn[2],]$Peakfile
        Pairname <- sprintf("%s:%s.vs.%s", Group, Donor1, Donor2)
        data.frame(Group, Donor1, Donor2, Sample1, Sample2, Peakfile1, Peakfile2, Pairname,
                   stringsAsFactors=FALSE)
    })[-1]
})
row.names(idr.jobs) <- idr.jobs$Pairname
idr.jobs <- within(idr.jobs, {
    jobname <- str_c("BC:", Pairname)
    logfile <- file.path(peakcall.dir, sprintf("%s.log", jobname))
    outbase <- file.path(idr.dir, make.names(Pairname))
    outfile <- str_c(outbase, "-npeaks-aboveIDR.txt")
    cmd <- mapply(makecmd, "time", "Rscript", basename(bca.script),
                  Peakfile1, Peakfile2, "-1", outbase, "0", "F", "p.value")
})

## Read all npeak tables
stopifnot(all(file.exists(idr.jobs$outfile)))
idr.npeaks <- ddply(idr.jobs, .(outfile), function(df) data.frame(df, setNames(read.table(df$outfile), c("S1", "S2", "IDR", "Npeak"))))
idr.npeaks <- data.frame(setNames(data.frame(do.call(rbind,str_split(idr.npeaks$Group, ":"))), c("Sampletype", "Celltype", "Timepoint")), idr.npeaks)
idr.npeaks <- within(idr.npeaks, {
    Celltype <- factor(Celltype, levels=c("Naive", "Memory"))
    Timepoint <- factor(Timepoint, levels=str_c("D", c("0", "1", "5", "14")))
    Comparison <- factor(str_extract(jobname, "\\d+\\.vs\\.\\d+$"))
    x <- do.call(rbind,str_split(Comparison, "\\.vs\\."))
    Donor1 <- factor(x[,1])
    Donor2 <- factor(x[,2])
    rm(x)
})

idr.npeak.summary <- ddply(idr.npeaks, .(Group, IDR), function(df) {
    data.frame(df[1,c("Sampletype", "Celltype", "Timepoint", "Group", "IDR")],
               rbind(summary(df$Npeak)))
})

saveRDS(idr.npeaks, file.path(idr.dir, "idr-npeaks.RDS"))
saveRDS(idr.npeak.summary, file.path(idr.dir, "idr-npeak-summary.RDS"))

idr.npeak.donor.mean <- ddply(idr.npeaks, .(Group, IDR), function(df) {
    alldonors <- sort(unique(c(levels(df$Donor1), levels(df$Donor2))))
    vsDonor <- llply(alldonors, function(d) df$Donor1 == d | df$Donor2 == d)
    donorMean <- laply(vsDonor, function(x) mean(df$Npeak[x]))
    otherMean <- laply(vsDonor, function(x) mean(df$Npeak[!x]))
    mr <- donorMean / otherMean
    x <- data.frame(df[1,c("Sampletype", "Celltype", "Timepoint", "Group", "IDR")],
                    Donor=alldonors, DonorMean=donorMean, OtherMean=otherMean, MeanRatio=mr)
    x[!is.na(x$MeanRatio),]
})

pdf("IDR plots.pdf", width=10, height=10)
ggplot(idr.npeak.summary) +
    geom_line(aes(x=IDR, y=Mean, ymax=Max., ymin=Min.), data=idr.npeak.summary) +
    geom_ribbon(aes(x=IDR, y=Mean, ymax=Max., ymin=Min.), data=idr.npeak.summary,
                alpha=0.1) +
    facet_grid(Sampletype + Celltype ~ Timepoint) +
    geom_line(aes(x=IDR, y=Npeak, group=jobname, colour=Comparison),
              data=idr.npeaks) +
    ggtitle("Npeak threshold vs IDR cutoff") + xlab("IDR threshold") + ylab("Number of peaks passing threshold")
dev.off()

pdf("IDR minmax plot.pdf")
ggplot(idr.npeak.summary) +
    aes(x=IDR, y=log2(Max./Min.), group=Group, colour=Sampletype) +
    geom_hline(linetype=2, size=2, yintercept=1) + geom_line() + ylim(0,3) +
    ggtitle("Npeak Max/Min ratio vs IDR threshold")
dev.off()

pdf("IDR meanratio plot.pdf", width=10, height=10)
ggplot(idr.npeak.donor.mean) +
    aes(x=IDR, y=log2(MeanRatio), group=Donor, colour=Donor) +
    geom_line() +
    facet_grid(Sampletype + Celltype ~ Timepoint) +
    ggtitle("Mean ratio of Npeaks involving Donor to other Npeaks") + xlab("IDR threshold") + ylab("Mean Npeak ratio")
dev.off()

pdf("IDR mean Npeak vs time.pdf", width=10, height=10)
ggplot(idr.npeak.summary) +
    aes(x=Timepoint, y=Mean, colour=factor(IDR), group=factor(IDR)) +
    geom_line() + scale_colour_discrete(name="IDR") +
    facet_grid(Sampletype ~ Celltype) +
    ggtitle("Mean Npeak vs Timepoint") + ylab("Mean Npeaks across all pairwise comparisons")
dev.off()

pdf("IDR median Npeak vs time.pdf", width=10, height=10)
ggplot(idr.npeak.summary) +
    aes(x=Timepoint, y=Median, colour=factor(IDR), group=factor(IDR)) +
    geom_line() + scale_colour_discrete(name="IDR") +
    facet_grid(Sampletype ~ Celltype) +
    ggtitle("Median Npeak vs Timepoint") + ylab("Median Npeaks across all pairwise comparisons")
dev.off()

pdf("IDR max Npeak vs time.pdf", width=10, height=10)
ggplot(idr.npeak.summary) +
    aes(x=Timepoint, y=Max., colour=factor(IDR), group=factor(IDR)) +
    geom_line() + scale_colour_discrete(name="IDR") +
    facet_grid(Sampletype ~ Celltype) +
    ggtitle("Max Npeak vs Timepoint") + ylab("Max Npeaks across all pairwise comparisons")
dev.off()

for (idr.threshold in sort(unique(idr.npeak.summary$IDR))) {
    idr.npeak.thresholds <- idr.npeak.summary[idr.npeak.summary$IDR == idr.threshold, c("Group", "Max.")]
    names(idr.npeak.thresholds)[2] <- "Npeaks"
    idr.npeak.thresholds$Peakfile <- file.path(peakcall.dir, sprintf("%s_peaks_sorted.encodePeak", make.names(idr.npeak.thresholds$Group)))
    idr.npeak.thresholds$Filtered.Peakfile <- file.path(peakcall.dir, sprintf("%s_peaks_IDR_%.2f_filtered.encodePeak", make.names(idr.npeak.thresholds$Group), idr.threshold))
    stopifnot(all(file.exists(idr.npeak.thresholds$Peakfile)))

    ## Produce filtered output files
    llply(seq(nrow(idr.npeak.thresholds)), function(i) {
        attach(idr.npeak.thresholds[i,])
        on.exit(detach())
        x <- read.narrowPeak(Peakfile, nrow=Npeaks)
        x$name <- basename(x$name)
        write.narrowPeak(x, Filtered.Peakfile)
        i
    }, .parallel=FALSE)
}

EOF <- NULL
EOF
