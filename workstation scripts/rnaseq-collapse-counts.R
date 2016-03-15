#!/usr/bin/Rscript

source("rnaseq-common2.R", chdir=TRUE)

library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

tsmsg("Reading raw counts")
sexp <- readRDS("rnaseq-counts-raw.RDS")
expdata <- colData(sexp)

tsmsg("Collapsing technical replicates of the same biological condition")
collapse.counts <- function(counts, by) {
    by <- as.factor(by)
    counts.by.biological.rep <- llply(levels(by), function(lev) {
        counts[, by == lev, drop=FALSE]
    })
    collapsed.counts <- llply(counts.by.biological.rep, rowSums)
    newcounts.matrix <- do.call(cbind, collapsed.counts)
    colnames(newcounts.matrix) <- levels(by)
    newcounts.matrix
}
collapsed.counts <- collapse.counts(assays(sexp)$counts, by=expdata$Sampletype)
collapsed.expdata <- expdata[!duplicated(expdata$Sampletype), c("Sampletype", "Donor", "Celltype", "Timepoint")]
rownames(collapsed.expdata) <- collapsed.expdata$Sampletype
## IMPORTANT: Reorder expdata to match column order of collapsed counts
collapsed.expdata <- collapsed.expdata[colnames(collapsed.counts),]
collapsed.sexp <- SummarizedExperiment(assays=list(counts=collapsed.counts),
                               rowData=rowData(sexp),
                               colData=as(collapsed.expdata, "DataFrame"))

tsmsg("Computing conditional quantile normalization offsets")
gn.len <- get.gene.lengths(TxDb.Hsapiens.UCSC.hg19.knownGene)
gn.gc <- get.gene.gc(TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Hsapiens.UCSC.hg19)
cqn.res <- cqn(assays(collapsed.sexp)$counts, x=gn.gc, lengths=gn.len)
assays(collapsed.sexp)$cqn.offset <- cqn.res$offset

tsmsg("Saving count data")
saveRDS(collapsed.sexp, "rnaseq-counts.RDS")
write.table(assays(collapsed.sexp)$counts, file="rnaseq-counts.txt")
tsmsg("Counting complete")
