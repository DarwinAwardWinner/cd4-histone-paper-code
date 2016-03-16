#!/bin/sh
# -*- mode:R -*-
#PBS
#PBS -j oe -o peak-count.log
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#SBATCH -c8
module load R/3.1.0;
cd /gpfs/home/rcthomps/Projects/sarah-cd4

Rscript --version
Rscript - <<'EOF';

setwd("/gpfs/home/rcthomps/Projects/sarah-cd4")

library(GenomicRanges)
library(Rsamtools)
library(Rsubread)
library(doParallel)
library(parallel)
library(plyr)
library(stringr)
library(xlsx)
registerDoParallel(cores=8)
## Not sure which of these is the real option. Some code seems to
## use different ones.
options(cores=8)
options(mc.cores=8)
library(BiocParallel)

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

gr.to.saf <- function(gr) {
    data.frame(Chr=as.vector(seqnames(gr)),
               Start=start(gr),
               End=end(gr),
               Strand=as.vector(strand(gr)),
               GeneID=names(gr))
}

## Calls featureCounts with appropriate options for peak counting
do.peakcount <- function(bam, saf, nthreads=getOption("cores"), ...) {
    if (is(saf, "GRanges"))
        saf <- gr.to.saf(saf)
    x <- featureCounts(
        bam, annot.ext=saf, isPairedEnd=FALSE,
        read2pos=5, strandSpecific=0,
        nthreads=nthreads, ...)$counts
    if (rownames(x)[1] == "geneid") {
        x <- x[-1,]
    }
    x
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

tsmsg("Reading sample data")
expdata <- fixdfchar(read.xlsx("sampledata.xlsx", 1))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
rownames(expdata) <- expdata$Sample
expdata$BAM.abspath <- file.path("/gpfs/home/rcthomps/Projects/sarah-cd4/bam_files", expdata$BAM)
stopifnot(all(file.exists(expdata$BAM.abspath)))

## Read merged peaks
peakcall.dir <- "/gpfs/home/rcthomps/Projects/sarah-cd4/peak_calls"
marks <- setdiff(levels(expdata$Sampletype), "input")
peak.files <- setNames(file.path(peakcall.dir, str_c(marks, "_peaks_IDR_filtered.encodePeak")),
                       marks)
stopifnot(all(file.exists(peak.files)))

incl.metacols <- c("Sample", "Donor", "Celltype", "Sampletype", "Timepoint", "Run", "Lane")

peaks <- llply(peak.files, read.narrowPeak, .parallel=TRUE)

peak.counts <- llply(marks, function(mark) {
  tsmsg("Counting reads in peaks for histone mark: ", mark)
  include.samples <- expdata$Sampletype %in% c(mark, "input")
  bamfiles <- setNames(expdata$BAM.abspath, make.names(rownames(expdata), unique=TRUE))[include.samples]

  tsmsg("Counting unique reads in peaks")
  counts <- do.peakcount(bamfiles, peaks[[mark]])
  sexp <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               colData=as(droplevels(expdata[include.samples,incl.metacols]), "DataFrame"),
                               rowData=peaks[[mark]])
  sexp
})

names(peak.counts) <- marks

tsmsg("Saving count data")
saveRDS(peak.counts, "peak-counts.RDS")
tsmsg("Counting complete")

## Terminate the heredoc, while remaining valid R code
EOF <- NULL
EOF
