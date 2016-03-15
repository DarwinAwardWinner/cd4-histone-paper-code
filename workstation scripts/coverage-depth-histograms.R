#!/usr/bin/env Rscript

library(dplyr)
library(magrittr)
library(stringr)
library(Biobase)
library(RColorBrewer)
library(rtracklayer)
library(GenomicRanges)

get.coverage.depth.counts <- function(BWFile) {
    x <- import(BWFile)
    stopifnot(is(x, "GRanges"))
    stopifnot(is.numeric(x$score))
    df <- data_frame(Width=as.numeric(width(x)), Depth=x$score)
    depth.factor <- df$Depth %>% .[.>0] %>% {1/min(.)}
    df$IntegerDepth <- round(df$Depth * depth.factor)
    y <- df %>% group_by(IntegerDepth) %>%
        do(data_frame(Count=as.numeric(sum(.$Width)))) %>%
        ungroup
    y$IntegerDepth %>% setdiff(0:max(.), .) %>%
        data_frame(IntegerDepth=., Count=0) %>%
        rbind(y, .) %>% arrange(IntegerDepth) %>%
        mutate(Depth=IntegerDepth / depth.factor,
               CumulativeCount=Count %>% cumsum)
}

sample.table <-
    list.files("bigwig_files/sample_cpm", pattern=".*\\.bw$", full.names = TRUE) %>%
    str_match("bigwig_files/sample_cpm/(((.*)\\.(.*)\\.(.*))\\.(.*))_normalized_cpm_pileup\\.bw$") %>%
    plyr::alply(2) %>% as_data_frame %>%
    set_colnames(c("BWFile", "Sample", "Group", "Sampletype", "Celltype", "Timepoint", "Donor")) %>%
    select(-BWFile, BWFile)

covcounts <- sample.table %>% group_by(BWFile) %>% do(merge(., get.coverage.depth.counts(.$BWFile)))

{
    plots <- covcounts %>% group_by(Sampletype) %>% do(loess100={
        ggplot(.) +
            aes(x=Depth, y=(Count+1) / sum(Count) * 100, group=Sample,
                color=Group, fill=Group,
                linetype=Donor, shape=Donor) +
            geom_smooth(method="loess", span=0.05) +
            scale_y_log10(name="Percent of Genome Covered at Depth") + xlim(0,100) +
            xlab("Depth per million reads") +
            ggtitle(sprintf("Count of Bases by Coverage Depth for %s\n(loess smoothed)", .$Sampletype[1]))
    }, loess25={
        ggplot(.) +
            aes(x=Depth, y=(Count+1) / sum(Count) * 100, group=Sample,
                color=Group, fill=Group,
                linetype=Donor, shape=Donor) +
            geom_smooth(method="loess", span=0.05) +
            scale_y_log10(name="Percent of Genome Covered at Depth") + xlim(0,25) +
            xlab("Depth per million reads") +
            ggtitle(sprintf("Count of Bases by Coverage Depth for %s\n(loess smoothed)", .$Sampletype[1]))
    }, line25= {
        ggplot(.) +
            aes(x=Depth, y=(Count+1) / sum(Count) * 100, group=Sample,
                color=Group, fill=Group,
                linetype=Donor, shape=Donor) +
            geom_line() +
            scale_y_log10(name="Percent of Genome Covered at Depth") + xlim(0,25) +
            xlab("Depth per million reads") +
            ggtitle(sprintf("Count of Bases by Coverage Depth for %s\n(exact count)", .$Sampletype[1]))
    })
    pdf("results/coverage-depth-histograms.pdf")
    print(plots$loess100)
    print(plots$loess25)
    print(plots$line25)
    dev.off()
}
