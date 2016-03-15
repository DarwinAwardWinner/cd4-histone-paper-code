#!/usr/bin/Rscript

source("common.R")
source("rnaseq-common.R")
source("promoter-data-common.R")

rownames(sampledata) <- with(sampledata, sprintf("%s-%s-%s-%s", donor, celltype, timepoint, sampletype))
sampledata$interaction <- with(sampledata, celltype:timepoint)
sampledata$interaction.withstype <- with(sampledata, sampletype:interaction)

## All together
run.aqm(sampledata, sample.counts,
        feature.subset=rowMax(sample.counts) > 5,
        reporttitle="Array Quality Metrics Report for Sarah CD4 ChIP-Seq Data",
        outdir=file.path("aqm", "chipseq-all-with-input"),
        intgroup="interaction.withstype",
        force=TRUE)

## All together minus input
run.aqm(sampledata, sample.counts,
        feature.subset=rowMax(sample.counts) > 5,
        reporttitle="Array Quality Metrics Report for Sarah CD4 ChIP-Seq Data",
        outdir=file.path("aqm", "chipseq-all"),
        intgroup="interaction.withstype",
        sample.subset=sampledata$sampletype != "input",
        force=TRUE)

## Separate for each mark
run.aqm(sampledata, sample.counts,
        feature.subset=rowMax(sample.counts) > 5,
        reporttitle.template="Array Quality Metrics Report for Sarah CD4 ChIP-Seq %s Data",
        outdir.template="aqm/chipseq-%s",
        intgroup="interaction",
        splitcols="sampletype",
        force=TRUE)
