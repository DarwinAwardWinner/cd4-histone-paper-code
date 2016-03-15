#!/usr/bin/Rscript

source("promoter-data-common2.R")
source("rnaseq-common2.R", chdir=TRUE)
sampledata <- readRDS("sampledata.RDS")
sample.counts <- readRDS("sample_counts.RDS")

selected.samples <- sampledata$Sampletype %in% c("input", "H3K4me2", "H3K4me3")
sampledata <- droplevels(sampledata[selected.samples,])
sample.counts <- sample.counts[,selected.samples]

make.all.group.contrasts <- function(groupnames, ref=1) {
    stopifnot(length(groupnames) >= 2)
    stopifnot(length(ref) == 1 && ref %in% seq_along(groupnames))
    sprintf("%s - %s", groupnames[-ref], groupnames[ref])
}

addsummary <- function(tables) {
    summary <- data.frame(diffcount=sapply(tables, nrow))
    c(list(summary=summary), tables)
}

design <- make.group.design(group=with(sampledata, Sampletype:Celltype:Timepoint), block=sampledata["Donor"])

input.groups <- colnames(design)[str_detect(colnames(design), "^input")]
all.timepoints <- levels(sampledata$Timepoint)
nonzero.timepoints <- setdiff(all.timepoints, "T0")

input.consistency.tests <- list(input.consistency=make.all.group.contrasts(input.groups))

timepoint.anova.tests <- foreach(stype=c("H3K4me2", "H3K4me3"), .combine=c) %:% foreach(ctype=c("Naive", "Memory"), .combine=c) %do% {
    x <- sprintf("%s.%s.%s - %s.%s.T0", stype, ctype, nonzero.timepoints, stype, ctype)
    names(x) <- sprintf("%s.%s.T0v%s", stype, ctype, nonzero.timepoints)
    setNames(nm=sprintf("%s.%s.AllT", stype, ctype), list(x))
}
timepoint.individual.tests <- unlist(unname(timepoint.anova.tests))

celltype.allt.tests <- foreach(stype=c("H3K4me2", "H3K4me3"), .combine=c) %do% {
    x <- sprintf("%s.Memory.%s - %s.Naive.%s", stype, all.timepoints, stype, all.timepoints)
    names(x) <- sprintf("%s.NvM.%s", stype, all.timepoints)
    setNames(nm=sprintf("%s.NvM.AllT", stype), list(x))
}
celltype.singlet.tests <- unlist(unname(celltype.allt.tests))

factorial.allt.tests <- foreach(stype=c("H3K4me2", "H3K4me3"), .combine=c) %do% {
    x <- sprintf("(%s.Memory.%s - %s.Naive.%s) - (%s.Memory.T0 - %s.Naive.T0)", stype, nonzero.timepoints, stype, nonzero.timepoints, stype, stype)
    names(x) <- sprintf("%s.Fac.%s", stype, nonzero.timepoints)
    setNames(nm=sprintf("%s.Fac.AllT", stype), list(x))
}
factorial.singlet.tests <- unlist(unname(factorial.allt.tests))

mi.vs.nf.tests <-
    setNames(sprintf("%s.Memory.T0 - %s.Naive.T336",
                     c("H3K4me2", "H3K4me3"), c("H3K4me2", "H3K4me3")),
             sprintf("%s.M0vN336",
                     c("H3K4me2", "H3K4me3")))

alltests <- c(input.consistency.tests,
              timepoint.anova.tests,
              timepoint.individual.tests,
              celltype.allt.tests,
              celltype.singlet.tests,
              factorial.allt.tests,
              factorial.singlet.tests,
              mi.vs.nf.tests)
