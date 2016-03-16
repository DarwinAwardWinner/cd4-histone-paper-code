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

options(error=function() traceback(2))

setwd("/gpfs/home/rcthomps/Projects/sarah-cd4")
source("chipseq-counting-common.R")

tsmsg("Reading sample data")
expdata <- fixdfchar(read.xlsx("sampledata.xlsx", 1))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
rownames(expdata) <- expdata$Sample
expdata$BAM.abspath <- file.path("/gpfs/home/rcthomps/Projects/sarah-cd4/bam_files", expdata$BAM)
stopifnot(all(file.exists(expdata$BAM.abspath)))

## Read merged peaks
peakcall.dir <- "/gpfs/home/rcthomps/Projects/sarah-cd4/peak_calls"
marks <- setdiff(levels(expdata$Sampletype), "input")
peak.files <- setNames(file.path(peakcall.dir, str_c(marks, "_peaks_IDR_0.25_filtered.encodePeak")),
                       marks)
stopifnot(all(file.exists(peak.files)))

incl.metacols <- c("Sample", "Donor", "Celltype", "Sampletype", "Timepoint", "Run", "Lane")

peaks <- llply(peak.files, read.narrowPeak, .parallel=TRUE)

peak.counts <- llply(marks, function(mark) {
  tsmsg("Counting reads in peaks for histone mark: ", mark)
  include.samples <- expdata$Sampletype %in% c(mark, "input")
  bamfiles <- setNames(expdata$BAM.abspath, rownames(expdata))[include.samples]

  tsmsg("Counting unique reads in peaks")
  counts <- countFragmentOverlapsByBedtools(bamfiles, peaks[[mark]], fraglength=1)[,]
  sexp <- SummarizedExperiment(assays=SimpleList(counts=counts),
                               colData=as(expdata[include.samples,incl.metacols], "DataFrame"),
                               rowData=peaks[[mark]])

  ## sexp <- summarizeOverlaps(features=peaks[[mark]], reads=bamfiles,
  ##                           mode=FivePrimeEnd, ignore.strand=TRUE,
  ##                           inter.feature=FALSE)
  colData(sexp) <- as(expdata[include.samples,incl.metacols], "DataFrame")
  sexp
})
names(peak.counts) <- marks

tsmsg("Saving count data")
saveRDS(peak.counts, "peak-counts.RDS")
tsmsg("Counting complete")

## Terminate the heredoc, while remaining valid R code
EOF <- NULL
EOF
