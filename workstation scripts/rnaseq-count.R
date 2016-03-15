#!/usr/bin/Rscript

source("rnaseq-common.R", chdir=TRUE)

library(annotate)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

tsmsg("Reading experiment data")
expdata <- read.xlsx("expdata_rnaseq.xlsx", 1)
expdata$Sample.ID <- as.character(expdata$Sample.ID)
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
expdata$Timepoint <- factor(expdata$Timepoint, levels=c("T0", "T24", "T120", "T336"))
# expdata$Sample.ID <- with(expdata, sprintf("R47L%s_E%s_%s.1mm", Lane, Barcode, Barcode_seq))
rownames(expdata) <- expdata$Sample.ID
expdata$Bam.File <- sprintf("tophat_results/%s/accepted_hits.bam", expdata$Sample.ID)
# expdata$Counts.File <- sprintf("htseq_counts/%s/counts.txt", expdata$Sample.ID)

bamfiles <- BamFileList(setNames(expdata$Bam.File, rownames(expdata)),
                        yieldSize=100000)
## Add expdata as metadata on the bam files
mcols(bamfiles) <- as(expdata, "DataFrame")

tsmsg("Reading annotation data")
gene.exons <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
## Add gene annotations as metadata on the exons-by-gene list
gene.annot <-
  DataFrame(row.names=names(gene.exons),
            ENTREZ=names(gene.exons),
            llply(c("SYMBOL", "GENENAME", "UNIGENE",
                    "ENSEMBL", "REFSEQ", "UCSCKG", "UNIPROT"),
                  function(i) {
                    CharacterList(lookUp(names(gene.exons),
                                         data="org.Hs.eg",
                                         what=i))
                  }))

mcols(gene.exons) <- gene.annot

tsmsg("Counting unique reads in genes")
counts <- summarizeOverlaps(gene.exons, bamfiles,
                            mode=IntersectionNotEmpty, ignore.strand=TRUE)
colData(counts) <- expdata

## tsmsg("Collapsing technical replicates")
## collapse.counts <- function(counts, by) {
##     newcounts <- do.call(cbind, llply(levels(by), function(i) rowSums(counts[,by==i, drop=FALSE])))
##     colnames(newcounts) <- levels(by)
##     newcounts
## }
## collapsed.counts <- collapse.counts(assays(counts)$counts, by=expdata$Sampletype)
## collapsed.expdata <- expdata[!duplicated(expdata$Sampletype), c("Sampletype", "Donor", "Celltype", "Timepoint")]
## rownames(collapsed.expdata) <- collapsed.expdata$Sampletype
## counts <- SummarizedExperiment(assays=list(counts=collapsed.counts),
##                                rowData=rowData(counts),
##                                colData=as(collapsed.expdata, "DataFrame"))

## tsmsg("Computing conditional quantile normalization offsets")
## gn.len <- get.gene.lengths(TxDb.Hsapiens.UCSC.hg19.knownGene)
## gn.gc <- get.gene.gc(TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Hsapiens.UCSC.hg19)
## cqn.res <- cqn(assays(counts)$counts, x=gn.gc, lengths=gn.len)
## assays(counts)$cqn.offset <- cqn.res$offset

tsmsg("Saving count data")
saveRDS(counts, "rnaseq-counts-raw.RDS")
write.table(assays(counts)$counts, file="rnaseq-counts-raw.txt")
tsmsg("Counting complete")
