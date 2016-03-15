#!/usr/bin/Rscript

source("rnaseq-diff-prep.R")

## Read batch-corrected norm factors
nf <- readRDS("batchnorm.RDS")

edger.analysis <-
    do.edger.analysis(counts, design, alltests,
                      feature.annot = gene.annot,
                      feature.subset = rowMax(counts) > 5,
                      group = with(expdata, Celltype:Timepoint),
                      group.blockfactors = expdata["Donor"],
                      sample.normfactors = nf,
                      .parallel = TRUE)

all.topresults.ql <- llply(edger.analysis$topresults.ql, function(x) {
    reorder.columns(x, start=c("geneSymbol", "description"))
})

all.topresults.lr <- llply(edger.analysis$topresults.lr, function(x) {
    reorder.columns(x, start=c("geneSymbol", "description"))
})

diffcounts.ql <- data.frame(cbind(diffcount=sapply(all.topresults.ql, nrow)))
diffcounts.lr <- data.frame(cbind(diffcount=sapply(all.topresults.lr, nrow)))

tsmsg("Saving results")
save.image("rnaseq-edger-results3.RDa")
write.xlsx.multisheet(addsummary(all.topresults.ql), "results/rnaseq-edger-topgenes3-ql.xlsx")
write.xlsx.multisheet(addsummary(all.topresults.lr), "results/rnaseq-edger-topgenes3-lr.xlsx")
tsmsg("Done.")
