#!/bin/sh
# -*- mode:R -*-
#PBS
#PBS -j oe -o promoter-group-count-2.5kb.log
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#SBATCH -c8
module load R/3.1.0;
cd /gpfs/home/rcthomps/Projects/sarah-cd4

Rscript --version
Rscript - <<'EOF';

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
    x <- featureCounts(bam, annot.ext=saf, isPairedEnd=FALSE, ...)$counts
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
## Get transcripts, annotated with ENTREZ gene ID
gene.tx <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
tx <- unlist(gene.tx)
tx$ENTREZ <- rep(factor(names(gene.tx)), elementLengths(gene.tx))
## Get 2500bp flanking regions around TSS
promoter.radius <- 2500
promoter <- flank(tx, start=TRUE, width=promoter.radius, both=TRUE)
## Group by overlap
promoter.reduced <- reduce(promoter)
promoter$overlap.group <- factor(findOverlaps(promoter, promoter.reduced, select="arbitrary"))
promoter$group <- factor(str_c(promoter$ENTREZ, promoter$overlap.group, sep=":"))
promoter.groups <- split(promoter, promoter$group)
## Add UCSC, ENTREZ IDs, and ENTREZ-derived gene-level annotations to groups
mcols(promoter.groups)$UCSCKG <- CharacterList(split(promoter$tx_name, promoter$group))
mcols(promoter.groups)$ENTREZ <- promoter$ENTREZ[cumsum(elementLengths(promoter.groups))]
## Delete ENTREZ from inner GRanges to avoid conflict
mcols(promoter.groups@unlistData)$ENTREZ <- NULL
group.annot <-
  DataFrame(row.names=names(promoter.groups),
            mcols(promoter.groups)[c("UCSCKG", "ENTREZ")],
            llply(c("SYMBOL", "GENENAME", "UNIGENE",
                    "ENSEMBL", "REFSEQ", "UNIPROT"),
                  function(i) {
                    CharacterList(lookUp(as.character(mcols(promoter.groups)$ENTREZ),
                                         data="org.Hs.eg",
                                         what=i))
                  }, .parallel=TRUE))
mcols(promoter.groups) <- group.annot

tsmsg("Counting unique reads in gene promoters")

counts <- do.count(setNames(expdata$BAM.abspath, rownames(expdata)), grl.to.saf(reduce(promoter.groups)),
                   nthreads=8, strandSpecific=0)
colnames(counts) <- rownames(expdata)

tsmsg("Constructing SummarizedExperiment object")
expdata <- expdata[c("Sample", "Donor", "Celltype", "Sampletype", "Timepoint", "Run", "Lane")]

sexp <- SummarizedExperiment(assays=SimpleList(counts=counts),
                             colData=as(expdata, "DataFrame"),
                             rowData=promoter.groups)
dimnames(sexp) <- dimnames(counts)

tsmsg("Saving count data")
saveRDS(sexp, "promoter-group-counts-2.5kb.RDS")
write.table(assays(sexp)$counts, file="promoter-group-counts-2.5kb.txt")
tsmsg("Counting complete")

## Terminate the heredoc, while remaining valid R
EOF <- NULL
EOF
