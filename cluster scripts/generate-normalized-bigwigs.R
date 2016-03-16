#!/bin/sh
# -*- mode:R -*-
#PBS
#PBS -j oe -o bigbin-count.log
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#SBATCH -c8
module load R/3.1.1;
cd /gpfs/home/rcthomps/Projects/sarah-cd4

Rscript --version

Rscript - <<'EOF';

source("chipseq-counting-common.R")
library(rtracklayer)
library(edgeR)
library(stringr)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)
library(lazyeval)
library(BiocParallel)

bprowapply <- function(DF, EXPR, BPPARAM=bpparam(), ...) {
    lazy.expr <- lazy(EXPR)
    args <- c(as.list(DF),
              list(FUN=function(...) lazyeval::lazy_eval(lazy.expr, data=list(...)),
                   BPPARAM=BPPARAM,
                   ...))
    do.call(bpmapply, args)
}

get.bam.coverage <- function(bam, fraglength=NULL) {
    aln <- readGAlignments(bam) %>% as("GRanges")
    if (!is.null(fraglength)) {
        aln <- resize(aln, width=fraglength)
    }
    aln %>% coverage
}

importBigWigAsRleList <- function(bwfile) {
    bwfile %>% import(format="BIGWIG") %>% coverage(weight=.$score)
}

## Attempt to delete file on failure
exportAtomic <- function(object, con, format, ...) {
    success <- FALSE
    if (is.character(con)){
        on.exit(if(!success) {
                    file.remove(con)
                })
    }
    result <- export(object, con, format, ...)
    success <- TRUE
    result
}

full.sampledata <- fixdfchar(read.xlsx("sampledata.xlsx", 1))
full.sampledata$Celltype <- factor(full.sampledata$Celltype, levels=c("Naive", "Memory"))
rownames(full.sampledata) <- full.sampledata$Sample
full.sampledata$BAM.abspath <- file.path("/gpfs/home/rcthomps/Projects/sarah-cd4/bam_files", full.sampledata$BAM)
stopifnot(all(file.exists(full.sampledata$BAM.abspath)))

marks <- setdiff(levels(full.sampledata$Sampletype), "input")

incl.metacols <- c("Sample", "Donor", "Celltype", "Sampletype", "Timepoint", "Run", "Lane")

bigbin.sexp <- readRDS("bigbin-counts.RDS")
colData(bigbin.sexp) <- DataFrame(full.sampledata[colnames(bigbin.sexp),])
input.samples <- bigbin.sexp %>% colData %$% Sampletype %>% equals("input")
input.sexp <- bigbin.sexp[,input.samples]
chip.sexp <- bigbin.sexp[,!input.samples]
dge <- chip.sexp %>% assay %>% DGEList %>% calcNormFactors(doWeighting=FALSE)
sampledata <- colData(chip.sexp) %>% as.data.frame %>% subset(Sampletype != "input") %>%
    mutate(
        lib.size=dge$samples$lib.size,
        norm.factor=dge$samples$norm.factors,
        BigWigCountFile=file.path("bigwig_files", "sample_counts",
        sprintf("%s.%s.%s.%s_count_pileup.bw",
                Sampletype, Celltype, Timepoint, Donor)),
        BigWigCPMFile=file.path("bigwig_files", "sample_cpm",
        sprintf("%s.%s.%s.%s_normalized_cpm_pileup.bw",
                Sampletype, Celltype, Timepoint, Donor)),
        Group=sprintf("%s.%s.%s", Sampletype, Celltype, Timepoint))
sampledata %$% c(BigWigCountFile, BigWigCPMFile) %>% dirname %>% unique %>%
    sapply(dir.create, showWarnings=FALSE, recursive=TRUE)

sampledata %>% cbind(SampleNum=seq(nrow(.))) %>% bprowapply(system.time({
    tsmsg <- function(...) {
        message(date(), ": ", ...)
    }
    if (file.exists(BigWigCountFile) && file.exists(BigWigCPMFile)) {
        tsmsg("Already processed sample ", SampleNum)
    } else {
        library(rtracklayer)
        library(magrittr)
        library(GenomicRanges)
        library(GenomicAlignments)
        get.bam.coverage <- function(bam, fraglength=NULL) {
            aln <- readGAlignments(bam) %>% as("GRanges")
            if (!is.null(fraglength)) {
                aln <- resize(aln, width=fraglength)
            }
            aln %>% coverage
        }
        ## Attempt to delete file on failure
        exportAtomic <- function(object, con, format, ...) {
            success <- FALSE
            if (is.character(con)){
                on.exit({
                    if(!success) {
                        tsmsg("Deleting incomplete file ", con)
                        file.remove(con)
                    }
                })
            }
            result <- export(object, con, format, ...)
            success <- TRUE
            invisible(result)
        }

        tsmsg("Working on sample ", SampleNum)
        tsmsg("Reading ", BAM.abspath)
        count.cov <- get.bam.coverage(BAM.abspath, fraglength=147)
        tsmsg("Normalizing coverage")
        norm.lib.size <- lib.size * norm.factor
        cpm.cov <- count.cov * 1e6 / norm.lib.size
        tsmsg("Exporting ", BigWigCountFile)
        exportAtomic(count.cov, BigWigCountFile, format="BIGWIG")
        tsmsg("Exporting ", BigWigCPMFile)
        exportAtomic(cpm.cov, BigWigCPMFile, format="BIGWIG")
        tsmsg("Finished with sample ", SampleNum)
    }
}), BPPARAM=SerialParam())

## Compute mean group CPM BigWigs
sampledata %>% split(.$Group) %>% bplapply(. %$% {
    tsmsg <- function(...) {
        message(date(), ": ", ...)
    }
    Group <- Group[1]
    group.bwfile <- file.path("bigwig_files", "group_cpm",
                              sprintf("%s_normalized_cpm_pileup.bw", Group))
    if (file.exists(group.bwfile)) {
        tsmsg("Already processed group ", Group)
    } else {
        library(rtracklayer)
        library(magrittr)
        library(GenomicRanges)
        library(GenomicAlignments)
        importBigWigAsRleList <- function(bwfile) {
            bwfile %>% import(format="BIGWIG") %>% coverage(weight=.$score)
        }
        ## Attempt to delete file on failure
        exportAtomic <- function(object, con, format, ...) {
            success <- FALSE
            if (is.character(con)){
                on.exit({
                    if(!success) {
                        tsmsg("Deleting incomplete file ", con)
                        file.remove(con)
                    }
                })
            }
            result <- export(object, con, format, ...)
            success <- TRUE
            invisible(result)
        }

        tsmsg("Working on group ", Group)
        stopifnot(all(file.exists(BigWigCPMFile)))
        num.samples <- length(BigWigCPMFile)
        tsmsg("Reading coverage for ", num.samples, " samples")
        tsmsg("Reading ", BigWigCPMFile[1])
        total.cov <- importBigWigAsRleList(BigWigCPMFile[1])
        for (bwfile in BigWigCPMFile[-1]) {
            tsmsg("Reading ", bwfile)
            total.cov %<>% add(importBigWigAsRleList(bwfile))
        }
        tsmsg("Computing mean coverage")
        mean.cov <- total.cov / num.samples
        tsmsg("Exporting ", group.bwfile)
        dir.create(dirname(group.bwfile), FALSE, TRUE)
        exportAtomic(total.cov, group.bwfile, format="BIGWIG")
        tsmsg("Done with group ", Group)
    }
}, BPPARAM=BatchJobsParam())
