#!/usr/bin/env Rscript

library(reshape2)
library(magrittr)
library(multcomp) # For glht
library(RColorBrewer)
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

library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(annotate)
library(rtracklayer)

flatten.CharacterList <- function(x, sep=",") {
    tryCatch(x <- as(x, "CharacterList"),
             error=function(...) NULL)
    if (is(x, "CharacterList"))
        setNames(rtracklayer:::pasteCollapse(x, collapse=sep),
                 names(x))
    else
        x
}

fix.cl.in.df <- function(df) {
    df[] <- llply(df, flatten.CharacterList, .parallel=TRUE)
}

load.object.from.file <- function(datafile, varname) {
    within(list(), load(datafile))[[varname]]
}

geneLengths <- function(txdb, which="max") {
    which <- match.fun(which)
    exbt <- exonsBy(txdb, "tx")
    txbg <- transcriptsBy(txdb, "gene")
    gene.txids <- split(txbg@unlistData$tx_id, rep(names(txbg), elementLengths(txbg)))
    txlen <- sapply(width(exbt), sum)
    gene.txlens <- lapply(gene.txids, function(i) txlen[i])
    gene.txlen <- sapply(gene.txlens, which)
    gene.txlen
}

makebins <- function(x, binsize=250) {
    cut(rank(x), round(length(x) / binsize))
}

bin.vector <- function(x, binsize=250) {
    f <- makebins(x, binsize)
    split(x, f)
}

chips <- c("H3K4me2", "H3K4me3", "H3K27me3")
chip.radius <- c(H3K4me2=1000, H3K4me3=1000, H3K27me3=2500)

## Load CPGi data and compute overlaps
cpgi <- read.table("CpGi.Ext.bed", sep="\t")
cpgi <- GRanges(seqnames=cpgi[[1]], ranges=IRanges(start=cpgi[[2]]+1, end=cpgi[[3]]), strand="*")

tss <- promoters(transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene), upstream=0, downstream=1)
names(tss) <- tss$tx_name
tss$cgi.dist <- mcols(distanceToNearest(tss, cpgi, select="arbitrary"))$distance

promoter.peak.table <- read.table("~/Projects/sarah-cd4-ng/results/Promoter-Table-PeakDistance.csv", header=TRUE)
promoter.peak.table$ENTREZ <- factor(as.character(promoter.peak.table$ENTREZ))
promoter.peak.table$CGI.Dist <- tss[promoter.peak.table$tx_name]$cgi.dist

rnaseq.analysis <- load.object.from.file("rnaseq-edger-results3.RDa", "edger.analysis")
dge <- rnaseq.analysis$dge
genelen.in.kb <- geneLengths(TxDb.Hsapiens.UCSC.hg19.knownGene) / 1000
genes <- Reduce(intersect, list(rownames(dge), names(genelen.in.kb), promoter.peak.table$ENTREZ))
samplecpm <- cpm(dge)
samplefpkm <- samplecpm[genes,] / expandAsMatrix(genelen.in.kb[genes], dim=dim(samplecpm[genes,]))
groupcpm <- quantify.by.group.DGEList(dge, with(expdata, Celltype:Timepoint), block=expdata["Donor"])
groupcpm <- groupcpm[genes,]
genelen.in.kb <- genelen.in.kb[genes]
groupfpkm <- groupcpm[genes,] / expandAsMatrix(genelen.in.kb[genes], dim=dim(groupcpm[genes,]))
grouplogfpkm <- log2(groupfpkm)

## Write CPM and FPKM tables
write.xlsx.multisheet(list(FPKM=data.frame(groupfpkm),
                           CPM=data.frame(groupcpm)),
                      "results/RNASeq FPKM by group.xlsx", row.names=TRUE)
write.xlsx.multisheet(list(FPKM=data.frame(samplefpkm),
                           CPM=data.frame(samplecpm)),
                      "results/RNASeq FPKM by sample.xlsx", row.names=TRUE)


## Bin plots
binsize <- 250

bindata <- {
    foreach (ctype=levels(expdata$Celltype), .combine=rbind) %do% {
        foreach (tp=levels(expdata$Timepoint), .combine=rbind) %do% {
            rnagroup <- str_c(ctype, tp, sep=".")
            logfpkm <- grouplogfpkm[,rnagroup]
            binfac <- makebins(logfpkm)
            logfpkm.bins <- split(logfpkm, binfac)
            bin.mean.logfpkm <- sapply(logfpkm.bins, mean)
            foreach (chip=chips, .combine=rbind) %do% {
                chipgroup <- str_c(chip, ctype, tp, sep=".")
                gene.peak.dist <- sapply(split(promoter.peak.table[[chipgroup]], promoter.peak.table$ENTREZ), min)[genes]
                gene.peak.presence <- gene.peak.dist[genes] <= chip.radius[chip]
                pp.bins <- split(gene.peak.presence, binfac)
                bin.peak.frac <- sapply(pp.bins, mean)
                data.frame(AveLogFPKM=bin.mean.logfpkm, PeakFrac=bin.peak.frac,
                           ChIP=chip, Condition=rnagroup, ChIP.Condition=chipgroup)
            }
        }
    }
}

{
    pdf("~/Projects/sarah-cd4-ng/results/FPKM Bin Peak Plots.pdf")
    dlply(bindata, .(ChIP.Condition), function(df) {
        chip <- as.character(df$ChIP[1])
        rnagroup <- as.character(df$Condition[1])
        print(ggplot(df) +
                  aes(x=AveLogFPKM, y=PeakFrac) + geom_point() +
                      ggtitle(sprintf(
                          "Fractions of %s-Gene FPKM Bins with %s Peaks\nin Promoters (radius %s) in Condition %s",
                          binsize, chip, chip.radius[chip], rnagroup)) +
                              xlab("Mean log2(FPKM) of genes in Bin") + ylab("Fraction of Promoters with Peaks") +
                                  ylim(0,1))
    })
    dev.off()
}

write.xlsx.multisheet(dlply(bindata, .(ChIP.Condition), `[`, c("AveLogFPKM", "PeakFrac")),
                      "../sarah-cd4-ng/results/FPKM-Bin-Peak-data.xlsx", row.names=FALSE)

## Tests of correlation
chip.alternative <- c(H3K4me2="up", H3K4me3="up", H3K27me3="down")

alldata <- {
    foreach (ctype=levels(expdata$Celltype), .combine=rbind) %do% {
        foreach (tp=levels(expdata$Timepoint), .combine=rbind) %do% {
            rnagroup <- str_c(ctype, tp, sep=".")
            logfpkm <- grouplogfpkm[,rnagroup]
            foreach (chip=chips, .combine=rbind) %do% {
                chipgroup <- str_c(chip, ctype, tp, sep=".")
                gene.peak.dist <- sapply(split(promoter.peak.table[[chipgroup]], promoter.peak.table$ENTREZ), min)[genes]
                gene.peak.presence <- gene.peak.dist <= chip.radius[chip]
                gene.cgi.dist <- sapply(split(promoter.peak.table$CGI.Dist, promoter.peak.table$ENTREZ), min)[genes]
                gene.cgi.presence <- gene.cgi.dist <= chip.radius[chip]
                data.frame(ENTREZ=rownames(grouplogfpkm),
                           Peak=ifelse(gene.peak.presence, "Peak", "No Peak"),
                           CGI=ifelse(gene.cgi.presence, "CGI", "No CGI"),
                           logFPKM=logfpkm,
                           ChIP=chip,
                           Condition=rnagroup,
                           ChIP.Condition=chipgroup)
            }
        }
    }
}

pdf("~/Projects/sarah-cd4-ng/results/FPKM by Peak Boxplots.pdf",
    width=12, height=8)
print(ggplot(alldata) + facet_grid(ChIP ~ Condition) +
          aes(y=logFPKM, x=Peak) + geom_boxplot() +
              ylab("log2(FPKM)") +
                  ggtitle("FPKM Distribution Conditioned on Peak Presence"))
dev.off()

pdf("~/Projects/sarah-cd4-ng/results/FPKM by Peak Violin Plots.pdf",
    width=12, height=8)
print(ggplot(alldata) + facet_grid(ChIP ~ Condition) +
          aes(y=logFPKM, x=Peak) + geom_violin() +
              ylab("log2(FPKM)") +
                  ggtitle("FPKM Distribution Conditioned on Peak Presence"))
dev.off()

## QQ Plots
qqdata <- ddply(alldata, .(ChIP, Condition), function(df) {
    xy <- split(df$logFPKM, df$Peak)
    df <- data.frame(qqplot(x=xy$`No Peak`, y=xy$Peak, plot.it=FALSE))
    names(df) <- c("NoPeak", "Peak")
    df
})

pdf("~/Projects/sarah-cd4-ng/results/FPKM by Peak QQ Plots.pdf",
    width=12, height=6)
print(ggplot(qqdata) + facet_grid(ChIP ~ Condition) +
          aes(x=NoPeak, y=Peak) + geom_line() + coord_equal() +
              geom_abline(slope=1, intercept=0, linetype="dotted") +
                  ggtitle("log2(FPKM) QQ Plot, Peak vs No Peak"))
dev.off()

## Different tests for differences:
## t-test with unequal variance
## Kolmogorov–Smirnov test for equal distributions
## Kendall correlation test
## Mann-Whitney U test
pvals <- ddply(alldata, .(ChIP, Condition), function(df) {
    xy <- split(df$logFPKM, df$Peak)

    chip.alternative <- c(H3K4me2="greater", H3K4me3="greater", H3K27me3="less")
    t.tst <- t.test(xy$Peak, xy$`No Peak`, alternative=chip.alternative[df$ChIP[1]])
    ks.tst <- ks.test(xy$`No Peak`, xy$Peak, alternative=chip.alternative[df$ChIP[1]])
    kend.tst <- cor.test(~ logFPKM + as.numeric(Peak), method="kendall", continuity=TRUE, data=df)

    chip.alternative <- c(H3K4me2="up", H3K4me3="up", H3K27me3="down")
    mwu.p <- with(df, geneSetTest(Peak == "Peak", logFPKM, alternative = chip.alternative[ChIP[1]]))

    data.frame(AveLogFPKMNoPeak=t.tst$estimate[2],
               AveLogFPKMPeak=t.tst$estimate[1],
               MedLogFPKMNoPeak=median(xy$`No Peak`),
               MedLogFPKMPeak=median(xy$Peak),
               AveDiffFPKM=t.tst$estimate[1] - t.tst$estimate[2],
               ConfInt95=t.tst$conf.int[ifelse(t.tst$alternative == "greater", 1, 2)],
               Kendall.Tau=kend.tst$estimate["tau"],
               Tstat=t.tst$statistic,
               TT.PValue=t.tst$p.value,
               KS.PValue=ks.tst$p.value,
               Kendall.PValue=kend.tst$p.value,
               MWU.PValue=mwu.p)
}, .parallel=TRUE)

write.xlsx(pvals, "../sarah-cd4-ng/results/FPKM-Peak-dependence-tests.xlsx",
                 row.names=FALSE)

## Linear model with peak presence as predictors and logFPKM as response
lmdata <- alldata
lmdata$Peak <- lmdata$Peak == "Peak"
lmdata <- dcast(lmdata, ENTREZ + Condition + logFPKM ~ ChIP, value.var="Peak")
## lmdata$Peak.Status <- str_c(, sep=".")
lmdata$CGI <- sapply(split(alldata$CGI, alldata$ENTREZ), function(x) "CGI" %in% x)[as.character(lmdata$ENTREZ)]

all.pairwise.test <- function(x, f, testfun=t.test, ...) {
    groups <- split(x, f)
    group.pairs <- alply(combn(names(groups), 2), 2, identity) %>%
        setNames(sapply(., str_c, collapse="_vs_"))
    llply(group.pairs, function(pair) {
        testfun(groups[[pair[1]]], groups[[pair[2]]], ...)
    })
}

## H3K4me3 and H3K27me3
{
    peak.status <- do.call(str_c, c(mapply(function(chip, presence) ifelse(presence, chip, "X"),
                                           chips[-1], lmdata[chips[-1]], SIMPLIFY=FALSE),
                                    list(sep=".")))
    peak.status[peak.status == "X.X"] <- "None"
    peak.status <- str_replace_all(peak.status, "X\\.|\\.?X$", "")
    peak.status <- relevel(factor(peak.status), ref="None")
    lmdata$PeakStatus <- peak.status

    difftests <-
        list(`KS Tests`=daply(lmdata, .(Condition), . %$% all.pairwise.test(logFPKM, PeakStatus, ks.test) %>% sapply(. %$% p.value)),
             `t-tests`=daply(lmdata, .(Condition), . %$% all.pairwise.test(logFPKM, PeakStatus, t.test) %>% sapply(. %$% p.value)))
    difftests.fdr <- difftests %>% lapply(. %>% {.[,] <- p.adjust(., "BH"); .})

    write.xlsx.multisheet(c(PValue=difftests, FDR=difftests.fdr), "/home/ryan/Projects/sarah-cd4-ng/results/FPKM by Peak Status Difference Tests.xlsx")

    pdf("~/Projects/sarah-cd4-ng/results/FPKM by Peak Status Boxplots.pdf",
        width=12, height=8)
    print(ggplot(lmdata) + facet_wrap(~ Condition, nrow=2) +
              aes(y=logFPKM, x=PeakStatus) + geom_boxplot() +
                  ylab("log2(FPKM)") +
                      ggtitle("FPKM Distribution Conditioned on Peak Status") + coord_flip())
    ## Same as above, but without bivalent genes
    print(ggplot(lmdata %>% subset(PeakStatus != "H3K4me3.H3K27me3")) + facet_wrap(~ Condition, nrow=2) +
              aes(y=logFPKM, x=PeakStatus) + geom_boxplot() +
                  ylab("log2(FPKM)") +
                      ggtitle("FPKM Distribution Conditioned on Peak Status") + coord_flip())
    dev.off()

    pdf("~/Projects/sarah-cd4-ng/results/FPKM by Peak Status Violin Plots.pdf",
        width=12, height=8)
    print(ggplot(lmdata) + facet_wrap(~ Condition, nrow=2) +
              aes(y=logFPKM, x=PeakStatus) + geom_boxplot() + geom_violin(alpha=0.25) +
                  ylab("log2(FPKM)") +
                      ggtitle("FPKM Distribution Conditioned on Peak Status") + coord_flip())
    ## Same as above, but without bivalent genes
    print(ggplot(lmdata %>% subset(PeakStatus != "H3K4me3.H3K27me3")) + facet_wrap(~ Condition, nrow=2) +
              aes(y=logFPKM, x=PeakStatus) + geom_boxplot() + geom_violin(alpha=0.25) +
                  ylab("log2(FPKM)") +
                      ggtitle("FPKM Distribution Conditioned on Peak Status") + coord_flip())
    dev.off()

    fits <- dlply(lmdata, .(Condition), function(df) {
        lm(logFPKM ~ H3K4me3 * H3K27me3, data=df)
    }, .parallel=TRUE)
}

## H3K4me2 and H3K4me3
{
    peak.status <- do.call(str_c, c(mapply(function(chip, presence) ifelse(presence, chip, "X"),
                                           chips[-3], lmdata[chips[-3]], SIMPLIFY=FALSE),
                                    list(sep=".")))
    peak.status[peak.status == "X.X"] <- "None"
    peak.status <- str_replace_all(peak.status, "X\\.|\\.?X$", "")
    peak.status <- relevel(factor(peak.status), ref="None")
    lmdata$PeakStatus <- peak.status

    difftests <-
        list(`KS Tests`=daply(lmdata, .(Condition), . %$% all.pairwise.test(logFPKM, PeakStatus, ks.test) %>% sapply(. %$% p.value)),
             `t-tests`=daply(lmdata, .(Condition), . %$% all.pairwise.test(logFPKM, PeakStatus, t.test) %>% sapply(. %$% p.value)))
    difftests.fdr <- difftests %>% lapply(. %>% {.[,] <- p.adjust(., "BH"); .})
    write.xlsx.multisheet(c(PValue=difftests, FDR=difftests.fdr), "/home/ryan/Projects/sarah-cd4-ng/results/FPKM by Peak Status Difference Tests H3K4.xlsx")

    pdf("~/Projects/sarah-cd4-ng/results/FPKM by Peak Status Boxplots H3K4.pdf",
        width=12, height=8)
    print(ggplot(lmdata) + facet_wrap(~ Condition, nrow=2) +
              aes(y=logFPKM, x=PeakStatus) + geom_boxplot() +
                  ylab("log2(FPKM)") +
                      ggtitle("FPKM Distribution Conditioned on Peak Status") + coord_flip())
    dev.off()

    pdf("~/Projects/sarah-cd4-ng/results/FPKM by Peak Status Violin Plots H3K4.pdf",
        width=12, height=8)
    print(ggplot(lmdata) + facet_wrap(~ Condition, nrow=2) +
              aes(y=logFPKM, x=PeakStatus) + geom_boxplot() + geom_violin(alpha=0.25) +
                  ylab("log2(FPKM)") +
                      ggtitle("FPKM Distribution Conditioned on Peak Status") + coord_flip())
    dev.off()

    fits <- dlply(lmdata, .(Condition), function(df) {
        lm(logFPKM ~ H3K4me3 * H3K27me3, data=df)
    }, .parallel=TRUE)
}

compute.ma <- function(..., do.log=FALSE) {
    within(xy.coords(...),{
        if (do.log) {
            x <- log2(x)
            y <- log2(y)
        }
        A <- (x + y)/2
        M <- y - x
    })
}

get.mareg.data <- function(lmdata, cond1, cond2, chipcond=cond1) {
    xt <- subset(lmdata, Condition == cond1, select=c(ENTREZ, logFPKM))
    yt <- subset(lmdata, Condition == cond2, select=c(ENTREZ, logFPKM))
    chipt <- subset(lmdata, Condition == chipcond, select=c(ENTREZ, H3K4me2, H3K4me3, H3K27me3, CGI))
    chipt$PeakStatus <- lapply(c("H3K4me2", "H3K4me3", "H3K27me3"), . %>% ifelse(chipt[[.]], ., "")) %>%
        {do.call(str_c, c(., list(sep=".")))} %>%
            str_replace_all("\\.\\.*", ".") %>%
                str_replace_all("^\\.\\.*|\\.\\.*$", "") %>%
                    str_replace("^$", "None") %>%
                        factor
    fullt <- merge(xt, yt, by="ENTREZ")
    fullt <- merge(fullt, chipt, by="ENTREZ")
    data.frame(fullt, compute.ma(fullt[c("logFPKM.x", "logFPKM.y")])[c("M", "A")])
}

quartiles <- summary(lmdata$logFPKM)[c(2,3,5)]

naive.madata <- get.mareg.data(lmdata, "Naive.D0", "Naive.D1")
memory.madata <- get.mareg.data(lmdata, "Memory.D0", "Memory.D1")

plotMAReg <- function(madata, groups) {
    madata$Group <- groups
    madata <- madata[!is.na(madata$Group),]
    ggplot(madata) + aes(x=A, y=M) +
        stat_density2d(geom="tile", aes(fill = ..density..), n=512, contour = FALSE) +
            geom_point(size=0.5) +
            scale_fill_gradientn(colours=c("white", brewer.pal(9, "Blues"))) +
                geom_smooth(method="lm") +
                    xlim(-10, 10) + ylim(-10, 10) +
                        facet_wrap(~ Group)
}

{
    pdf("~/Projects/sarah-cd4-ng/results/MA Plots by Peak status.pdf", width=14, height=8)
    ## H3K4me3 & H3K27me3 vs neither
    with(naive.madata, {
        groups <- factor(rep(NA, nrow(naive.madata)),
                         levels=c("Neither H3K4me3 nor H3K27me3",
                         "Bivalent H3K4me3 and H3K27me3"))
        groups[!H3K4me3 & !H3K27me3] <- "Neither H3K4me3 nor H3K27me3"
        groups[H3K4me3 & H3K27me3] <- "Bivalent H3K4me3 and H3K27me3"
        print(plotMAReg(naive.madata, groups) +
              ggtitle("Bivalent at D0 vs expression MA plot at D1 in Naive"))
    })
    with(memory.madata, {
        groups <- factor(rep(NA, nrow(memory.madata)),
                         levels=c("Neither H3K4me3 nor H3K27me3",
                         "Bivalent H3K4me3 and H3K27me3"))
        groups[!H3K4me3 & !H3K27me3] <- "Neither H3K4me3 nor H3K27me3"
        groups[H3K4me3 & H3K27me3] <- "Bivalent H3K4me3 and H3K27me3"
        print(plotMAReg(memory.madata, groups) +
              ggtitle("Bivalent at D0 vs expression MA plot at D1 in Memory"))
    })
    ## H3K4me3 and no H3K27me3 vs neither
    with(naive.madata, {
        groups <- factor(rep(NA, nrow(naive.madata)),
                         levels=c("Neither H3K4me3 nor H3K27me3",
                         "H3K4me3 & No H3K27me3"))
        groups[!H3K4me3 & !H3K27me3] <- "Neither H3K4me3 nor H3K27me3"
        groups[H3K4me3 & !H3K27me3] <- "H3K4me3 & No H3K27me3"
        print(plotMAReg(naive.madata, groups) +
              ggtitle("H3K4me3 at D0 vs expression MA plot at D1 in Naive"))
    })
    with(memory.madata, {
        groups <- factor(rep(NA, nrow(naive.madata)),
                         levels=c("Neither H3K4me3 nor H3K27me3",
                         "H3K4me3 & No H3K27me3"))
        groups[!H3K4me3 & !H3K27me3] <- "Neither H3K4me3 nor H3K27me3"
        groups[H3K4me3 & !H3K27me3] <- "H3K4me3 & No H3K27me3"
        print(plotMAReg(memory.madata, groups) +
              ggtitle("H3K4me3 at D0 vs expression MA plot at D1 in Memory"))
    })
    ## H3K4me3 and no H3K27me3 vs neither
    with(naive.madata, {
        groups <- factor(rep(NA, nrow(naive.madata)),
                         levels=c("Neither H3K4me3 nor H3K27me3",
                         "H3K27me3 & No H3K4me3"))
        groups[!H3K4me3 & !H3K27me3] <- "Neither H3K4me3 nor H3K27me3"
        groups[!H3K4me3 & H3K27me3] <- "H3K27me3 & No H3K4me3"
        print(plotMAReg(naive.madata, groups) +
              ggtitle("H3K27me3 at D0 vs expression MA plot at D1 in Naive"))
    })
    with(memory.madata, {
        groups <- factor(rep(NA, nrow(naive.madata)),
                         levels=c("Neither H3K4me3 nor H3K27me3",
                         "H3K27me3 & No H3K4me3"))
        groups[!H3K4me3 & !H3K27me3] <- "Neither H3K4me3 nor H3K27me3"
        groups[!H3K4me3 & H3K27me3] <- "H3K27me3 & No H3K4me3"
        print(plotMAReg(memory.madata, groups) +
              ggtitle("H3K27me3 at D0 vs expression MA plot at D1 in Memory"))
    })
    dev.off()
}

naive.madata %<>% transform(ExpBin=cut(A, breaks=c(-Inf, -4, 1, Inf), labels=c("Low", "Med", "High")))
memory.madata %<>% transform(ExpBin=cut(A, breaks=c(-Inf, -4, 1, Inf), labels=c("Low", "Med", "High")))

## Do t-tests on Low, Med, and High expression
{
    naive.matests <- with(naive.madata, {
        ps <- ifelse(H3K4me3,
                     ifelse(H3K27me3, "Bivalent", "H3K4me3"),
                     ifelse(H3K27me3, "H3K27me3", "Neither")) %>% factor
        f <- ExpBin:ps
        levels(f) %<>% make.names
        design <- model.matrix(~0 + f)
        colnames(design) <- levels(f)
        form <- as.formula(str_c("M ~ 0 + ", str_c(colnames(design), collapse="+")))
        fit <- lm(form, data=data.frame(design))
        ct <- unlist(llply(levels(ExpBin), function(explev) {
            sprintf(c("%s.Bivalent - %s.Neither",
                      "%s.H3K4me3 - %s.Neither",
                      "%s.H3K27me3 - %s.Neither"),
                    explev, explev)
        }))
        cmat <- t(makeContrasts(contrasts=ct, levels=design))
        htests <- glht(fit, cmat)
        individual.tests <- summary(htests, adjusted("none"))
        df <- data.frame(individual.tests$test[c("coefficients", "sigma", "tstat", "pvalues")]) %>%
            setNames(c("Estimate", "Std. Error", "t value", "PValue"))
        rownames(df) %<>% str_c("Naive ", .)
        df
    })
    memory.matests <- with(memory.madata, {
        ps <- ifelse(H3K4me3,
                     ifelse(H3K27me3, "Bivalent", "H3K4me3"),
                     ifelse(H3K27me3, "H3K27me3", "Neither")) %>% factor
        f <- ExpBin:ps
        levels(f) %<>% make.names
        design <- model.matrix(~0 + f)
        colnames(design) <- levels(f)
        form <- as.formula(str_c("M ~ 0 + ", str_c(colnames(design), collapse="+")))
        fit <- lm(form, data=data.frame(design))
        ct <- unlist(llply(levels(ExpBin), function(explev) {
            sprintf(c("%s.Bivalent - %s.Neither",
                      "%s.H3K4me3 - %s.Neither",
                      "%s.H3K27me3 - %s.Neither"),
                    explev, explev)
        }))
        cmat <- t(makeContrasts(contrasts=ct, levels=design))
        htests <- glht(fit, cmat)
        individual.tests <- summary(htests, adjusted("none"))
        df <- data.frame(individual.tests$test[c("coefficients", "sigma", "tstat", "pvalues")]) %>%
            setNames(c("Estimate", "Std. Error", "t value", "PValue"))
        rownames(df) %<>% str_c("Memory ", .)
        df
    })
    matests <- rbind(naive.matests, memory.matests)
    matests$FDR <- p.adjust(matests$PValue, "BH")
    matests <- matests[as.vector(matrix(seq_along(rownames(matests)), byrow=TRUE, ncol=3)),]
    write.xlsx(matests, "../sarah-cd4-ng/results/MA-peak-tests.xlsx",
               row.names=TRUE)
}
