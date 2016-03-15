#!/usr/bin/Rscript

source("rnaseq-diff-prep.R")

expdata$Batch <- factor(ifelse(expdata$Timepoint %in% c("T0", "T120"), "B1", "B2"))

design <- model.matrix(~Batch, data=expdata)
design.intonly <- model.matrix(~1, data=expdata)

dimnames(counts) <- list(gene=sprintf("gene%s",seq(nrow(counts))),
                         sample=sprintf("sample%s", seq(ncol(counts))))
rownames(design) <- NULL
dge <- DGEList(counts)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design=design, prior.df=0)
save(list=c("dge", "design"), file="testdata.rda")

load("testdata.rda")
library(edgeR)

## batch.analysis <-
##     do.edger.analysis(counts, design,
##                       ## feature.subset = rowMax(counts) >= 50,
##                       disp.tagwise.opts = list(prior.df = 0),
##                       group=expdata$Batch)
## dge <- batch.analysis$dge

split.dge.into.quantiles <- function(dge, n=100) {
    logcpm <- dge$AveLogCPM
    qlogcpm <- quantile(unique(logcpm), seq(0,1,length.out=n+1))
    cutlogcpm <- cut(logcpm, qlogcpm, include.lowest=TRUE)
    ind <- split(seq(nrow(dge)), cutlogcpm)
    dgesplit <- lapply(ind, function(i) dge[i,])
}

nbin <- 100
genes.per.bin <- 20
dgequant <- split.dge.into.quantiles(dge, nbin)

smalldispgenes <-
    do.call(c, lapply(dgequant, function(d) {
        rownames(d)[order(d$tagwise.dispersion)[1:genes.per.bin]]
    }))

dgesmalldisp <- dge[rownames(dge) %in% smalldispgenes,]
dgesmalldisp <- calcNormFactors(dgesmalldisp)
dgesmalldisp <- estimateDisp(dgesmalldisp, design, prior.df = 0)
nf.smalldisp <- dgesmalldisp$samples$norm.factors

## edger.analysis <-
##     do.edger.analysis(counts, design,
##                       sample.normfactors = nf.smalldisp,
##                       feature.subset = rowQ(counts, 17) > 5,
##                       group=expdata$Batch,
##                       tests=list(Batch="BatchB2"))


std.nf <- dge$samples$norm.factors
smd.nf <- nf.smalldisp
nf.delta <- smd.nf - std.nf
qplot(x=std.nf, y=nf.delta, color=expdata$Batch)
qplot(x=std.nf, y=smd.nf, color=expdata$Batch) + geom_abline(slope=1, intercept=0)

saveRDS(smd.nf, "batchnorm.RDS")
