#!/usr/bin/Rscript

source("rnaseq-common2.R", chdir=TRUE)

addsummary <- function(tables) {
    summary <- data.frame(diffcount=sapply(tables, nrow))
    c(list(summary=summary), tables)
}

assignFrom(read.counts("rnaseq-counts.RDS"))
## Translate from "24 hours" to "1 day"
tp.nametrans <- c(T0="D0", T24="D1", T120="D5", T336="D14")
expdata$Timepoint <- factor(tp.nametrans[expdata$Timepoint], levels=tp.nametrans)

design <- make.group.design(with(expdata, Celltype:Timepoint), expdata["Donor"])

## Test definitions
celltypes <- levels(expdata$Celltype)
all.timepoints <- levels(expdata$Timepoint)
nonzero.timepoints <- setdiff(all.timepoints, "D0")

timepoint.anova.tests <- setNames(llply(celltypes, function(ct) {
    setNames(sprintf("%s.%s - %s.D0", ct, nonzero.timepoints, ct),
             sprintf("%s.D0v%s", ct, nonzero.timepoints))
}), nm=str_c(celltypes, ".AllT"))
timepoint.single.tests <- as.list(unlist2(timepoint.anova.tests))
celltype.singlet.tests <-
    as.list(setNames(sprintf("Memory.%s - Naive.%s", all.timepoints, all.timepoints),
                     sprintf("NvM.%s", all.timepoints)))
celltype.allt.test <- list(NvM.AllT=unlist(celltype.singlet.tests))
factorial.singlet.tests <-
    as.list(setNames(sprintf("(Memory.%s - Memory.D0) - (Naive.%s - Naive.D0)",
                             nonzero.timepoints, nonzero.timepoints),
                     sprintf("Fac.%s", nonzero.timepoints)))
factorial.allt.test <- list(Fac.AllT=unlist(factorial.singlet.tests))
mi.vs.nf.test <- list(MD0vND14="Memory.D0 - Naive.D14")

batch.tests <- list(B1vB2="((Naive.D1 + Memory.D1 + Naive.D14 + Memory.D14) - (Naive.D0 + Memory.D0 + Naive.D5 + Memory.D5)) / 4",
                    crossbatch="(Naive.D0 + Memory.D0 + Naive.D14 + Memory.D14 - (Naive.D1 + Memory.D1 + Naive.D5 + Memory.D5)) / 4")

donor.var.test <- list(InterDonor=sprintf("Donor%s", 1:3))
alltests <- c(timepoint.anova.tests, timepoint.single.tests,
              celltype.allt.test, celltype.singlet.tests,
              factorial.allt.test, factorial.singlet.tests,
              mi.vs.nf.test, batch.tests, donor.var.test)
