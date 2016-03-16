#!/bin/sh
# -*- mode:R -*-
#PBS
#PBS -j oe -o window-count.log
#PBS -l nodes=1:ppn=1,mem=15gb,walltime=48:00:00
#SBATCH -c8
module load R/3.1.0;
cd /gpfs/home/rcthomps/Projects/sarah-cd4

Rscript --version
Rscript - <<'EOF';

options(error=function() traceback(2))

setwd("/gpfs/home/rcthomps/Projects/sarah-cd4")
source("chipseq-counting-common.R")
library(Matrix)

tsmsg("Reading sample data")
expdata <- fixdfchar(read.xlsx("sampledata.xlsx", 1))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
rownames(expdata) <- expdata$Sample
expdata$BAM.abspath <- file.path("/gpfs/home/rcthomps/Projects/sarah-cd4/bam_files", expdata$BAM)
stopifnot(all(file.exists(expdata$BAM.abspath)))

marks <- setdiff(levels(expdata$Sampletype), "input")

incl.metacols <- c("Sample", "Donor", "Celltype", "Sampletype", "Timepoint", "Run", "Lane")

histone.width <- 148

tsmsg("Building sliding windows")
windows <- sliding.windows(BamFileList(expdata$BAM.abspath),
                           window.interval=50, window.size=histone.width)
print(length(windows))

window.counts <- {
    tsmsg("Counting reads in histone-sized windows")
    bamfiles <- setNames(expdata$BAM.abspath, make.unique(rownames(expdata)))
    ## counts <- featureCounts.fragments(bamfiles, windows, fraglength=histone.width,
    ##                                   bigmatrix.backing.file="window-counts.bigmat")
    counts <- countFragmentOverlapsByBedtools(bamfiles, windows, fraglength=histone.width,
                                              bigmatrix.backing.file="window-counts.bigmat")
    sexp <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                 colData=as(droplevels(expdata[,incl.metacols]), "DataFrame"),
                                 rowData=windows)
    sexp
}

tsmsg("Saving count data")
saveRDS(window.counts, "window-counts.RDS")
tsmsg("Counting complete")

## Terminate the heredoc, while remaining valid R code
EOF <- NULL
EOF
