#!/usr/bin/Rscript

library(stringr)
library(dplyr)
library(ggplot2)
library(Biobase)
library(RColorBrewer)
library(Matrix)
library(rtracklayer)
library(magrittr)

## idr.npeak.summary <- readRDS("idr-npeak-summary.RDS") %>%
##     filter(Sampletype %in% c("H3K4me2", "H3K4me3"))
## ggplot(idr.npeak.summary) +
##     aes(x=Timepoint, y=Max., colour=factor(IDR), group=factor(IDR)) +
##     geom_line() + scale_colour_discrete(name="IDR") +
##     facet_grid(Sampletype ~ Celltype) +
##     ggtitle("Max Npeak vs Timepoint") + ylab("Max Npeaks across all pairwise comparisons")

idr.npeak.summary <- readRDS("idr-npeak-summary.RDS") %>%
    filter(Sampletype %in% c("H3K4me2", "H3K4me3"),
           IDR == 0.02) %>%
    droplevels

## Need cairo_pdf for the "≤" symbol to work
cairo_pdf("results/1kb/IDR max Npeak vs time.pdf", width=10, height=10)
ggplot(idr.npeak.summary) +
    aes(x=Timepoint, y=Max., colour=Sampletype,
        group=Sampletype:Celltype, linetype=Celltype, shape=Celltype) +
    geom_point(size=4) + geom_line(alpha=0.5, size=1) +
    scale_colour_discrete(name="Sample Type") +
    scale_linetype_manual(name="Cell Type", values=c("dashed","solid")) +
    scale_shape_discrete(name="Cell Type") +
    ylim(c(0,max(idr.npeak.summary$Max.))) +
    ggtitle("Max Npeak across groups") + ylab("Number of peaks with IDR≤2% in best donor pair")
dev.off()
