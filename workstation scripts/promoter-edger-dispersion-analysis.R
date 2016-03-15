#!/usr/bin/Rscript

source("promoter-diff-prep.R")

singlegroups <- split(seq(nrow(sampledata)),
                      str_replace_all(with(sampledata, Sampletype:Celltype:Timepoint), ":", "."))
singledesigns <- llply(singlegroups, function(group) cbind(Intercept=rep(1, length(group))))

singlegroup.edger.analyses <- llply(names(singlegroups), function(i) {
    do.edger.analysis(sample.counts[,singlegroups[[i]]], singledesigns[[i]], tests=NULL,
                      disp.tagwise.opts=list(prior.df=0))
}, .parallel=TRUE)

groups <- do.call(c, llply(sampledata[c("Sampletype","Celltype","Timepoint")],
                           function(x) split(seq(nrow(sampledata)), x)))
designs <- llply(groups, function(x) drop.zero.columns(design[x,]))

group.edger.analyses <- llply(names(groups), function(i) {
    do.edger.analysis(sample.counts[,groups[[i]]], designs[[i]], tests=NULL,
                      disp.tagwise.opts=list(prior.df=0))
})

disp <- sapply(group.edger.analyses, function(x) x$dge.lr$common.dispersion)
logCPM <- data.frame(llply(group.edger.analyses, function(x) x$dge.lr$logCPM))

get.bcv.limits <- function(x) {
    disps <- c(x$common.dispersion, x$trended.dispersion, x$tagwise.dispersion)
    bcvs <- sqrt(disps)
    bcvs <- na.omit(bcvs)
    c(min(bcvs), max(bcvs))
}

get.abundance.limits <- function(x) {
    logCPM <- na.omit(x$logCPM)
    c(min(logCPM), max(logCPM))
}

group.ylims <- Reduce(function(x, y) c(min(c(x[1],y[1])), max(c(x[2], y[2]))),
                llply(group.edger.analyses, function(x) get.bcv.limits(x$dge.lr)))
group.xlims <- Reduce(function(x, y) c(min(c(x[1],y[1])), max(c(x[2], y[2]))),
                llply(group.edger.analyses, function(x) get.abundance.limits(x$dge.lr)))

pdf("groupdisps.pdf")
for (i in names(group.edger.analyses)) {
    plotBCV(group.edger.analyses[[i]]$dge.lr, main=i,
            xlim=group.xlims, ylim=group.ylims)
}
dev.off()

singlegroup.ylims <- Reduce(function(x, y) c(min(c(x[1],y[1])), max(c(x[2], y[2]))),
                llply(singlegroup.edger.analyses, function(x) get.bcv.limits(x$dge.lr)))
singlegroup.xlims <- Reduce(function(x, y) c(min(c(x[1],y[1])), max(c(x[2], y[2]))),
                llply(singlegroup.edger.analyses, function(x) get.abundance.limits(x$dge.lr)))

pdf("singlegroupdisps.pdf")
for (i in names(singlegroup.edger.analyses)) {
    plotBCV(singlegroup.edger.analyses[[i]]$dge.lr, main=i,
            xlim=singlegroup.xlims, ylim=singlegroup.ylims)
}
dev.off()

save.image("promoter-edger-dispersion-analysis.RDa")
