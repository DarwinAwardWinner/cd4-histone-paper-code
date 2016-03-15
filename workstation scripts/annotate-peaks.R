source("rnaseq-common2.R", chdir=TRUE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rtracklayer)
library(plyr)
library(parallel)
library(doParallel)
registerDoParallel(cores=8)
## Not sure which of these is the real option. Some code seems to
## use different ones.
options(cores=parallel:::detectCores())
options(mc.cores=parallel:::detectCores())
library(BiocParallel)
library(stringr)

## X and Y should be GRanges, IRanges, etc.
select.nearest <- function(x, y) {
    y[nearest(x,y)]
}

annotate.by.transcript <- function(peaks, tdb, peakloc=c("near", "midpoint", "far", "summit")) {
    all.transcripts <- transcripts(tdb)
    nearest.transcripts <- select.nearest(peaks, all.transcripts)
    transcript.distance <- distance(peaks, nearest.transcripts, ignore.strand=TRUE)
    in.transcript <- Rle(transcript.distance == 0)
    kgID <- elementMetadata(nearest.transcripts)$tx_name
    DataFrame(genic=in.transcript, transcript.distance=transcript.distance, kgID=kgID)
}

annotate.by.cds <- function(peaks, tdb, peakloc=c("near", "midpoint", "far", "summit")) {
    all.cds <- cds(tdb)
    nearest.cds <- select.nearest(peaks, all.cds)
    cds.distance <- distance(peaks, nearest.cds, ignore.strand=TRUE)
    coding <- Rle(cds.distance == 0)
    DataFrame(coding=coding)
}

annotate.by.promoter <- function(peaks, tdb, promoter.size=1000,
                                 peakloc=c("near", "midpoint", "far", "summit")) {
    peakloc <- match.arg(peakloc)
  all.transcripts <- transcripts(tdb)
  tss.loc <- ifelse(strand(all.transcripts) == "+",
                    start(all.transcripts),
                    end(all.transcripts))
  all.tss <- all.transcripts
  start(all.tss) <- end(all.tss) <- tss.loc
  nearest.tss <- select.nearest(peaks, all.tss)
  tss.distance <- distance(peaks, nearest.tss, ignore.strand=TRUE)
  in.promoter <- Rle(tss.distance <= promoter.size)
  DataFrame(in.promoter=in.promoter, tss.distance=tss.distance, tss.kgID=elementMetadata(nearest.tss)$tx_name)
}

annotate.by.exon <- function(peaks, tdb) {
  all.exons <- exons(tdb)
  nearest.exons <- select.nearest(peaks, all.exons)
  exon.distance <- distance(peaks, nearest.exons, ignore.strand=TRUE)
  exonic <- Rle(exon.distance == 0)
  DataFrame(exonic=exonic)
}

dnase.file <- "dnase-hs-clusters/hg19-DNaseI-HS-clusters.bed"
dnase.clusters <- import(dnase.file, format="BED", asRangedData=FALSE, genome="hg19")
annotate.by.dnase.hs <- function(peaks, dnase.clusters) {
  nearest.clusters <- select.nearest(peaks, dnase.clusters)
  cluster.distance <- distance(peaks, nearest.clusters, ignore.strand=TRUE)
  DataFrame(DNaseI.hs.distance=cluster.distance)
}

granges.to.dataframe <- function(gr, ignore.strand=FALSE) {
  df.columns <- list(row.names=names(gr),
                     chr=seqnames(gr),
                     start=start(gr),
                     end=end(gr))
  if (!ignore.strand) {
    df.columns <- c(df.columns, list(strand=strand(gr)))
  }
  df.columns <- c(df.columns, as.list(elementMetadata(gr)))
  do.call(DataFrame, df.columns)
}

cbind.unique.colnames <- function(...) {
  items <- list(...)
  if (length(items) <= 1) {
    do.call(cbind, items)
  } else {
    combined.item <- items[[1]]
    for (item in items[-1]) {
      combined.item <- cbind(combined.item, item[! colnames(item) %in% colnames(combined.item)])
    }
    combined.item[!duplicated(names(combined.item))]
  }
}

infile <- "merged_peaks.bed"
peaks <- import(infile, format="BED", asRangedData=FALSE, genome="hg19")
## Migrate the name column to the names
names(peaks) <- elementMetadata(peaks)$name
elementMetadata(peaks)$name <- NULL

tdb <- get.ucsc.transcript.db("hg19", "knownGene")

annot <- do.call(cbind,
                 llply(list(annotate.by.promoter,
                            annotate.by.transcript,
                            annotate.by.cds,
                            annotate.by.exon),
                       do.call,
                       list(peaks, tdb)))
annot <- cbind(annot, annotate.by.dnase.hs(peaks, dnase.clusters))
annot$category <- factor(NA, levels=c("promoter", "exon", "intron", "UTR", "intergenic"))
annot$category[as.vector(!annot$genic)] <- "intergenic"
annot$category[as.vector(annot$exonic)] <- "exon"
annot$category[as.vector(annot$coding & !annot$exonic)] <- "intron"
annot$category[as.vector(annot$genic & !annot$coding)] <- "UTR"
annot$category[as.vector(annot$in.promoter)] <- "promoter"

gene.kgxref <- get.kgxref(annot$kgID)
tss.kgxref <- get.kgxref(annot$tss.kgID)
## Make sure TSS-related columns are marked as such
names(tss.kgxref) <- sprintf("tss.%s", names(tss.kgxref))
annot <- cbind.unique.colnames(annot, gene.kgxref, tss.kgxref)
elementMetadata(peaks) <- cbind.unique.colnames(elementMetadata(peaks), annot)
peaks.table <- as.data.frame(granges.to.dataframe(peaks))
saveRDS(peaks, "annotated_peaks.RDS")
write.csv(peaks.table, "annotated_peaks.csv")

## Seems to be incapable of writing an xlsx file this large.
## write.xlsx2(peaks.table, "annotated_peaks.xlsx")
