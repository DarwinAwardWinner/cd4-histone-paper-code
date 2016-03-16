#!/bin/sh
# -*- mode:R -*-
#PBS
#PBS -j oe -o promoter-gene-count.log
#SBATCH -c8
module load R/3.0.1;
cd $PBS_O_WORKDIR

Rscript - <<EOF;

library(annotate)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(Rsamtools)
library(xlsx)
library(plyr)
options(cores=8)
options(mc.cores=8)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(Rsubread)
## Not sure which of these is the real option. Some code seems to
## use different ones.
options(cores=parallel:::detectCores())
options(mc.cores=parallel:::detectCores())
registerDoParallel()


tsmsg <- function(...) {
  message(date(), ": ", ...)
}

## If llply is called on an unnamed character vector, use the
## character vector itself for the names.
llply <- function(.data, ...) {
  if (is.character(.data) && is.null(names(.data))) {
    names(.data) <- .data
  }
  plyr::llply(.data, ...)
}

fixdfchar <- function(df) {
    for (i in seq_along(df)) {
        if (is.factor(df[[i]]) && !any(duplicated(df[[i]]))) {
            df[[i]] <- as.character(df[[i]])
        }
    }
    df
}

## This merges exons into genes (GRanges to GRangesList)
gr.to.grl <- function(gr, featureType="exon", attrType="gene_id") {
    gr <- gr[gr$type %in% featureType]
    split(gr, as.character(mcols(gr)[[attrType]]))
}

## This converts a GRangesList into the SAF ("Simplified annotation
## format")
grl.to.saf <- function(grl) {
    gr <- unlist(grl)
    data.frame(Chr=as.vector(seqnames(gr)),
               Start=start(gr),
               End=end(gr),
               Strand=as.vector(strand(gr)),
               GeneID=rep(names(grl), elementLengths(grl)))
}

do.count <- function(bam, saf, ...) {
    x <- featureCounts(bam, file.type="BAM", annot.ext=saf,
                       isPairedEnd=FALSE,
                       ...)$counts
    if (rownames(x)[1] == "geneid") {
        x <- x[-1,]
    }
    x
}

tsmsg("Reading sample data")
expdata <- fixdfchar(read.xlsx("sampledata.xlsx", 1))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
rownames(expdata) <- expdata$Sample
expdata$BAM.abspath <- file.path("/gpfs/home/rcthomps/Projects/sarah-cd4/bam_files", expdata$BAM)
stopifnot(all(file.exists(expdata$BAM.abspath)))

tsmsg("Reading annotation data")
gene.tx <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
## Add gene annotations as metadata on the exons-by-gene list
promoter.radius <- 1000
gene.promoter <- flank(gene.tx, start=TRUE, width=promoter.radius, both=TRUE)
gene.annot <-
  DataFrame(row.names=names(gene.promoter),
            ENTREZ=names(gene.promoter),
            llply(c("SYMBOL", "GENENAME", "UNIGENE",
                    "ENSEMBL", "REFSEQ", "UCSCKG", "UNIPROT"),
                  function(i) {
                    CharacterList(lookUp(names(gene.promoter),
                                         data="org.Hs.eg",
                                         what=i))
                  }, .parallel=TRUE))
mcols(gene.promoter) <- gene.annot

tsmsg("Counting unique reads in gene promoters")

counts <- do.count(setNames(expdata$BAM.abspath, rownames(expdata)), grl.to.saf(reduce(gene.promoter)),
                   nthreads=8, strandSpecific=0)
colnames(counts) <- rownames(expdata)

tsmsg("Constructing SummarizedExperiment object")
expdata <- expdata[c("Sample", "Donor", "Celltype", "Sampletype", "Timepoint", "Run", "Lane")]

sexp <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             colData=as(expdata, "DataFrame"),
                             rowData=gene.promoter)
dimnames(sexp) <- dimnames(counts)

tsmsg("Saving count data")
saveRDS(sexp, "promoter-gene-counts.RDS")
write.table(assays(sexp)$counts, file="promoter-gene-counts.txt")
tsmsg("Counting complete")

## Terminate the heredoc, while remaining valid R
EOF <- NULL
EOF
