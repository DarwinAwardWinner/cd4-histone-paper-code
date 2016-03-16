#!/bin/sh
# -*- mode:R -*-
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -j oe -o rnaseq-count-star.log

## This makes sure we're running R 3.0.1.
module load R/3.0.1;
## This reports the version at the beginning of the log, to verify
## that this shell trickery worked
Rscript --version
Rscript - <<'EOF';

setwd(file.path(Sys.getenv("HOME"), "Projects", "sarah-restim"))

library(xlsx)
library(Rsubread)
library(rtracklayer)
library(annotate)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(plyr)
library(doParallel)
registerDoParallel(cores=8)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

star.results.dir <- "star_results"
star.results.file <- "Aligned.out.sam"
star.results.bam.file <- "Aligned.out.bam"

message("Reading experiment metadata")

expdata <- read.xlsx("expdata.xlsx", sheetName="RNASeq")

for (i in names(expdata)) {
    if (is.factor(expdata[[i]]) && nlevels(expdata[[i]]) == nrow(expdata)) {
        expdata[[i]] <- as.character(expdata[[i]])
    }
}
rownames(expdata) <- expdata$Sample
expdata <- within(expdata, {
    star.outdir <- file.path(star.results.dir, Sample)
    star.outfile <- file.path(star.outdir, star.results.bam.file)
    outfile.exists <- file.exists(star.outfile)
})
stopifnot(all(expdata$outfile.exists))

tsmsg("Reading annotation data")

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


tsmsg("Reading annotation data")
gene.exons <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
## Add gene annotations as metadata on the exons-by-gene list
metacols <- c("SYMBOL", "GENENAME", "UNIGENE", "ENSEMBL", "REFSEQ", "UCSCKG", "UNIPROT")
gene.annot <-
  DataFrame(row.names=names(gene.exons),
            ENTREZ=names(gene.exons),
            llply(metacols,
                  function(i) {
                    CharacterList(lookUp(names(gene.exons),
                                         data="org.Hs.eg",
                                         what=i))
                  }))
names(gene.annot)[-1] <- metacols

mcols(gene.exons) <- gene.annot

saf <- grl.to.saf(gene.exons)

do.count <- function(bam, saf, ...) {
    featureCounts(bam, annot.ext=saf,
                  isPairedEnd=FALSE,
                  ...)$counts
}

tsmsg("Counting unique sense strand fragments in genes")
sense.counts <- do.count(expdata[,"star.outfile"], saf,
                         nthreads=8, strandSpecific=1)
tsmsg("Counting unique antisense strand fragments in genes")
antisense.counts <- do.count(expdata[,"star.outfile"], saf,
                             nthreads=8, strandSpecific=2)
totalcounts <- sense.counts + antisense.counts

save.image("rnaseq-count.rda")

tsmsg("Constructing SummarizedExperiment object")
incl.metacols <- c("Sample", "Donor", "Celltype", "Stim", "Restim", "Run", "Lane", "Barcode")
expdata <- as(droplevels(expdata[,incl.metacols]), "DataFrame")

sexp <- SummarizedExperiment(assays=
                             SimpleList(counts=totalcounts,
                                        sense.counts=sense.counts,
                                        antisense.counts=antisense.counts),
                             colData=DataFrame(expdata),
                             rowData=gene.exons)

tsmsg("Saving count data")
saveRDS(sexp, "rnaseq-counts-star.RDS")
write.table(assays(sexp)$counts, file="rnaseq-counts-star.txt")
tsmsg("Counting complete")

EOF <- NULL
EOF
