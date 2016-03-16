#!/bin/sh
# -*- mode:R -*-
#PBS
#PBS -j oe -o bigbin-count.log
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#SBATCH -c8
module load R/3.1.0;
cd /gpfs/home/rcthomps/Projects/sarah-cd4

Rscript --version
Rscript - <<'EOF';

options(error=function() traceback(2))

setwd("/gpfs/home/rcthomps/Projects/sarah-cd4")
source("chipseq-counting-common.R")

tsmsg("Reading sample data")
expdata <- fixdfchar(read.xlsx("sampledata.xlsx", 1))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
rownames(expdata) <- expdata$Sample
expdata$BAM.abspath <- file.path("/gpfs/home/rcthomps/Projects/sarah-cd4/bam_files", expdata$BAM)
stopifnot(all(file.exists(expdata$BAM.abspath)))

marks <- setdiff(levels(expdata$Sampletype), "input")

incl.metacols <- c("Sample", "Donor", "Celltype", "Sampletype", "Timepoint", "Run", "Lane")

tsmsg("Creating 10kb bins")
bigbins <- sliding.windows(BamFileList(expdata$BAM.abspath),
                           window.interval=10000)

bigbin.counts <- {
    tsmsg("Counting reads in bigbins")
    bamfiles <- setNames(expdata$BAM.abspath, make.names(rownames(expdata), unique=TRUE))

    counts <- featureCounts.fiveprime(bamfiles, bigbins)$counts[,]
    sexp <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                 colData=as(droplevels(expdata[,incl.metacols]), "DataFrame"),
                                 rowData=bigbins)
    sexp
}

tsmsg("Saving count data")
saveRDS(bigbin.counts, "bigbin-counts.RDS")
tsmsg("Counting complete")

## Terminate the heredoc, while remaining valid R code
EOF <- NULL
EOF
