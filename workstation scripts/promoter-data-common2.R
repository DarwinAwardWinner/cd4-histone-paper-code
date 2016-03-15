#!/usr/bin/Rscript

source("common.R", chdir=TRUE)

library(rtracklayer)
library(DiffBind)
library(plyr)
library(stringr)
library(RColorBrewer)
library(xlsx)

tsmsg <- function(...) {
  message(date(), ": ", ...)
}

tsmsgf <- function(...) {
  tsmsg(sprintf(...))
}

retrieveParallelResult <- function(p) {
  result <- tryCatch(collect(p)[[1]], error=function(...) NULL)
  if (is(result, "try-error")) {
    stop(attr(result, "condition")$message)
  } else {
    if (is.null(result)) {
      invisible(result)
    } else {
      result
    }
  }
}

inSubprocess <- function(expr) {
  retrieveParallelResult(mcparallel(expr))
}

`%=%` <- function(lhs, rhs) {
  lhs <- as.character(as.symbol(substitute(lhs)))
  rhs.parallel <- mcparallel(rhs)
  delayedAssign(lhs, retrieveParallelResult(rhs.parallel), assign.env=parent.frame(1))
}

## Additional args are passed to every call of addDataFrame
write.xlsx.multisheet <- function(data.frames, file, sheetNames=names(data.frames), ...) {
  if (is.null(sheetNames)) {
    sheetNames <- str_c("Sheet", seq_along(data.frames))
  }
  ## Ensure correct number of sheetNames
  stopifnot(length(sheetNames) == length(data.frames))
  ## Fill in missing names if needed
  sheetNames[is.na(sheetNames)] <- str_c("Sheet", seq_along(data.frames))[is.na(sheetNames)]
  wb <- createWorkbook()
  sheets <- llply(sheetNames, function(x) createSheet(wb, sheetName=x))
  mlply(cbind(sheet=sheets, x=data.frames), .fun=addDataFrame, .parallel=FALSE, ...)
  saveWorkbook(wb, file)
}

renice.self <- function(niceness=19)  {
  system2("renice", args=c("-n", as.character(niceness), Sys.getpid()))
}

makeNiceCluster <- function(...) {
  cl <- makeCluster(...)
  clusterCall(cl, renice.self)
  cl
}

parse.fraglength.from.xls <- function(filename) {
  starting.text <- sapply(filename,
                          function (x) str_c(readLines(x, n=100),
                                             collapse="\n"))
  as.integer(str_match(starting.text, "\n# d = ([0-9]*)\n")[,2])
}

read.macs.xls <- function(...) {
  read.table(sep="\t", header=TRUE, comment.char="#", ...)
}

write.macs.xls <- function(...) {
  write.table(sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, ...)
}

## FDR is given in percent!
filter.macs.xls.by.fdr <- function(infile, outfile, fdr=0.1, fdr.percent=fdr*100) {
  x <- read.macs.xls(infile, check.names=FALSE)
  x <- x[x$`FDR(%)` <= fdr.percent,]
  write.macs.xls(x, outfile)
  x
}

select.nearest <- function(x, y) {
  y[nearest(x,y)]
}

countOverlapsFromBam <- function(query, bamfile, ...) {
  ## param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE), which=query,
  ##                       simpleCigar = TRUE, reverseComplement = TRUE, what = ShortRead:::.readAligned_bamWhat())
  ## subject <- as(readAligned(bamfile, type="BAM", param=param), "GRanges")
  inSubprocess(countOverlaps(query, readAlignedRanges(bamfile), ...))
}

## Read minimal information from BAM file to create a GRanges
readGRangesFromBam <- function(bamfile) {
  x <- scanBam(bamfile, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                          what=c("rname", "pos", "strand", "qwidth")))[[1]]
  GRanges(seqnames=x$rname,
          ranges=IRanges(start=x$pos, width=x$qwidth),
          strand=Rle(x$strand))
}

readMCSFromBam <- function(bamfile) {
  x <- scanBam(bamfile, param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                          what=c("rname", "pos", "strand")))[[1]]
  data.frame(chr=x$rname, pos=x$pos, strand=x$strand)
}

granges.to.mcs <- function(gr) {
  data.frame(chr=seqnames(gr),
             pos=ifelse(strand(gr) == "+", start(gr), end(gr)),
             strand=strand(gr))
}

## Modify the 3prime end of each range to be newwidth away from the 5prime end
adjust.width <- function(gr, newwidth) {
  pstrand <- strand(gr) == "+"
  ## Adjust end of plus strand tags
  end(gr[pstrand]) <- start(gr[pstrand]) + newwidth - 1
  ## Adjust start of minus strand tags
  start(gr[!pstrand]) <- end(gr[!pstrand]) - newwidth + 1
  gr
}

estimateDisp <- function(dge, design, common.disp=TRUE, trended.disp=TRUE, tagwise.disp=TRUE) {
  dge <- calcNormFactors(dge)
  if (common.disp)
    dge <- estimateGLMCommonDisp(dge, design, verbose=TRUE)
  if (trended.disp)
    dge <- estimateGLMTrendedDisp(dge, design)
  if (tagwise.disp)
    dge <- estimateGLMTagwiseDisp(dge, design)
  dge
}

expdata <- read.xlsx("expdata.xls", 1)
expdata$celltype <- factor(expdata$celltype, levels=c("N", "M"))
expdata
expdata$donor <- factor(as.character(expdata$donor))
expdata$ChIP.BAM <- file.path("aligned_reads", sprintf("%s.bam", expdata$chip))
expdata$Control.BAM <- file.path("aligned_reads", sprintf("%s.bam", expdata$input))
## Not real xls files
expdata$Peaks.Xls <- file.path("macs_output", expdata$expname, sprintf("%s_peaks.xls", expdata$expname))
expdata$Peaks.Xls.filtered <- file.path("macs_output", expdata$expname, sprintf("%s_peaks_0.1FDR.xls", expdata$expname))

## Redundant with above "xls" files
expdata$Peaks.Bed <- file.path("macs_output", expdata$expname, sprintf("%s_peaks.bed", expdata$expname))
expdata$fraglength <- parse.fraglength.from.xls(expdata$Peaks.Xls)

mean.fraglength <- mean(expdata$fraglength)

## sampledata.chip <- rename(cbind(expdata[c("chip", "donor", "celltype", "timepoint", "ChIP.BAM")],
##                                 sampletype=factor(sprintf("IP.%s", expdata$ab), levels=c("input", sprintf("IP.%s", levels(expdata$ab))))),
##                           c(chip="sample", ChIP.BAM="BAM"))
## sampledata.input <- rename(cbind(expdata[c("input", "donor", "celltype", "timepoint", "Control.BAM")],
##                                  sampletype=factor("input", levels=c("input", sprintf("IP.%s", levels(expdata$ab))))),
##                            c(input="sample", Control.BAM="BAM"))
## sampledata.input <- unique(sampledata.input)
## sampledata <- rbind(sampledata.chip, sampledata.input)
## rownames(sampledata) <- sampledata$sample

sampledata <- read.xlsx("sampledata.xlsx", 1)
rownames(sampledata) <- sampledata$Sample
sampledata$BAM <- as.character(file.path("aligned_reads", sprintf("%s.bam", sampledata$Sample)))
sampledata$Celltype <- factor(sampledata$Celltype, levels=c("Naive", "Memory"))
sampledata$Sampletype <- relevel(sampledata$Sampletype, ref="input")
sampledata$Timepoint <- factor(sampledata$Timepoint, levels=c("T0", "T24", "T120", "T336"))

saveRDS(sampledata, "sampledata.RDS")

## library(annotate)
## library(BSgenome.Hsapiens.UCSC.hg19)
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## library(org.Hs.eg.db)

## txbygene <- unlist(transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene"))
## gene.annot <-
##   DataFrame(row.names=names(txbygene),
##             ENTREZ=names(txbygene),
##             llply(c("SYMBOL", "GENENAME", "UNIGENE",
##                     "ENSEMBL", "REFSEQ", "UCSCKG", "UNIPROT"),
##                   function(i) {
##                     CharacterList(lookUp(names(txbygene),
##                                          data="org.Hs.eg",
##                                          what=i))
##                   }))

## Get promoter regions
tdb <- get.ucsc.transcript.db("hg19", "knownGene")
all.transcripts <- transcripts(tdb)
tss.loc <-
  ifelse(strand(all.transcripts) == "+",
         start(all.transcripts),
         end(all.transcripts))
## Eliminate duplicate TSS sites
tssdup <- duplicated(sprintf("%s:%s:%s",
                             seqnames(all.transcripts),
                             strand(all.transcripts),
                             as.vector(tss.loc)))
all.transcripts <- all.transcripts[!tssdup]
tss.loc <- as.vector(tss.loc[!tssdup])

all.promoters <- all.transcripts
## Start 1000 bp before, end 1000 bp after, but constrain to ends of
## chromosomes
start(all.promoters) <- pmax(tss.loc - 1000, 1)
end(all.promoters) <- pmin(tss.loc + 1000, seqlengths(all.promoters)[as.vector(seqnames(all.promoters))])

x <- as.vector(elementMetadata(all.promoters)$tx_name)
kgxref <- get.ucsc.table("knownGene","kgXref")
row.names(kgxref) <- make.unique(as.character(kgxref$kgID), sep="/")
promoter.annotations <-
  data.frame(row.names=x, kgxref[x,],
             TSS.coord=sprintf("%s:%s", seqnames(all.promoters), tss.loc))

## Read data
sample.counts <-
  tryCatch({
    sample.counts <- readRDS("sample_counts.RDS")
  }, error = function(e) {
    sample.counts <- llply(sampledata$BAM,
                           function(BAM)
                           countOverlaps(adjust.width(readGRangesFromBam(BAM), newwidth=mean.fraglength),
                                         query=all.promoters),
                           .parallel=TRUE)
    sample.counts <- do.call(cbind, sample.counts)
    colnames(sample.counts) <- rownames(sampledata)
    rownames(sample.counts) <- promoter.annotations$kgID
    saveRDS(sample.counts, "sample_counts.RDS")
    sample.counts
  })
