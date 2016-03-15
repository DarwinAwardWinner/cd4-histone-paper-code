#!/usr/bin/env/Rscript
source("rnaseq-common2.R", chdir=TRUE)

library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(annotate)
library(rtracklayer)
library(magrittr)
library(scales)

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

wanted.cols <- c("ENTREZ", "CpGi", "PValue.H3K4me2", "FDR.H3K4me2",
                 "logFC.H3K4me2", "PValue.H3K4me3", "FDR.H3K4me3", "logFC.H3K4me3",
                 "PValue.H3K27me3", "FDR.H3K27me3", "logFC.H3K27me3", "PValue.RNASeq",
                 "FDR.RNASeq", "logFC.RNASeq")

full.tabs <- readRDS("edger-merged-tables3-ql.RDS")
full.tabs <- llply(full.tabs, function(df) {
    df <- df[names(df) %in% wanted.cols]
    df$ENTREZ <- as.character(df$ENTREZ)
    df
})
## Load logcpm table
rnatabs <- load.object.from.file("rnaseq-edger-results3.RDa", "edger.analysis")$results.ql
logcpm <- rnatabs[[1]]
logcpm <- logcpm[str_detect(names(logcpm), "logCPM\\.")]
logcpm[] <- lapply(logcpm, function(x) ifelse(is.na(x), -Inf, x))
## names(logcpm) %<>% str_replace("logCPM\\.", "")
## Bin into low/med/hi expression with thresholds of 0 and 4
exprbins <- logcpm
exprbins[] <- lapply(exprbins, cut, breaks=c(-Inf, 0,4, Inf), labels=c("Low", "Med", "High"))
names(exprbins) %<>% str_replace("logCPM\\.", "exprbin.")

## save(list=c("wanted.cols", "full.tabs", "rnatabs", "logcpm", "exprbins"), file="~/Dropbox/temp.rda")

pdf("~/Dropbox/salomon/sarah-results/densityplot.pdf", width=12, height=8)
## X-axis: H3K4me3 N0 vs N5
## Y-axis: H3K27me3 N0 vs N1
## Bin by expr, plot density in separate panels for each bin
## OR: Plot points, color by expr
## Ditto for memory
{
    df <- data.frame(full.tabs$Naive.D0vD1["ENTREZ"],
                     logFC.H3K27me3.D1=full.tabs$Naive.D0vD1$logFC.H3K27me3,
                     logFC.H3K4me3.D5=full.tabs$Naive.D0vD5$logFC.H3K4me3)
    df$logCPM <- logcpm[df$ENTREZ,]$logCPM.Naive.D0
    df$logCPM[is.na(df$logCPM)] <- -Inf
    df$bin <- exprbins[df$ENTREZ,]$exprbin.Naive.D0
    df$bin[is.na(df$bin)] <- "Low"
    print(
        ggplot(df) +
            aes(x=logFC.H3K4me3.D5, y=logFC.H3K27me3.D1) +
            stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
            scale_fill_gradient(low="white", high=muted("blue")) +
            stat_density2d(colour="black", size=.25) + xlim(-2,2) + ylim(-2,2) +
            facet_wrap(~bin) + coord_equal() +
            ggtitle("Histone changes binned by Day 0 RNA expression, Naive")
    )

    df <- data.frame(full.tabs$Memory.D0vD1["ENTREZ"],
                     logFC.H3K27me3.D1=full.tabs$Memory.D0vD1$logFC.H3K27me3,
                     logFC.H3K4me3.D5=full.tabs$Memory.D0vD5$logFC.H3K4me3)
    df$logCPM <- logcpm[df$ENTREZ,]$logCPM.Memory.D0
    df$logCPM[is.na(df$logCPM)] <- -Inf
    df$bin <- exprbins[df$ENTREZ,]$exprbin.Memory.D0
    df$bin[is.na(df$bin)] <- "Low"
    print(
        ggplot(df) +
            aes(x=logFC.H3K4me3.D5, y=logFC.H3K27me3.D1) +
            stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
            scale_fill_gradient(low="white", high=muted("blue")) +
            stat_density2d(colour="black", size=.25) + xlim(-2,2) + ylim(-2,2) +
            facet_wrap(~bin) + coord_equal() +
            ggtitle("Histone changes binned by Day 0 RNA expression, Memory")
    )
}
dev.off()
system("pdfcrop ~/Dropbox/salomon/sarah-results/densityplot.pdf")
system("mv ~/Dropbox/salomon/sarah-results/densityplot-crop.pdf ~/Dropbox/salomon/sarah-results/densityplot.pdf")
