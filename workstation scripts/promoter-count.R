#!/usr/bin/Rscript
#SBATCH -c8
source("rnaseq-common2.R", chdir=TRUE)

library(annotate)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

tsmsg("Reading experiment data")
expdata <- read.xlsx("expdata_rnaseq.xlsx", 1)
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
expdata$Sample.ID <- with(expdata, sprintf("R47L%s_E%s_%s.1mm", Lane, Barcode, Barcode_seq))
rownames(expdata) <- expdata$Sample.ID
expdata$Bam.File <- sprintf("tophat_results/%s/accepted_hits.bam", expdata$Sample.ID)
expdata$Counts.File <- sprintf("htseq_counts/%s/counts.txt", expdata$Sample.ID)

bamfiles <- BamFileList(setNames(expdata$Bam.File, rownames(expdata)),
                        yieldSize=100000)
## Add expdata as metadata on the bam files
mcols(bamfiles) <- as(expdata, "DataFrame")

tsmsg("Reading annotation data")
gene.tx <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
## Add gene annotations as metadata on the exons-by-gene list
promoter.radius <- 1000
gene.promoter <- flank(gene.tx, width=promoter.radius, both=TRUE)
gene.annot <- 
  DataFrame(row.names=names(gene.promoter),
            ENTREZ=names(gene.promoter),
            llply(c("SYMBOL", "GENENAME", "UNIGENE",
                    "ENSEMBL", "REFSEQ", "UCSCKG", "UNIPROT"),
                  function(i) {
                    CharacterList(lookUp(names(gene.promoter),
                                         data="org.Hs.eg",
                                         what=i))
                  }))
mcols(gene.promoter) <- gene.annot

tsmsg("Counting unique reads in gene promoters")
counts <- summarizeOverlaps(gene.promoter, bamfiles,
                            mode=IntersectionNotEmpty, ignore.strand=TRUE)
colData(counts) <- as(expdata, "DataFrame")

tsmsg("Computing conditional quantile normalization offsets")
gn.len <- get.gene.lengths(TxDb.Hsapiens.UCSC.hg19.knownGene)
gn.gc <- get.gene.gc(TxDb.Hsapiens.UCSC.hg19.knownGene, BSgenome.Hsapiens.UCSC.hg19)
cqn.res <- cqn(assays(counts)$counts, x=gn.gc, lengths=gn.len)
assays(counts)$cqn.offset <- cqn.res$offset

tsmsg("Saving count data")
saveRDS(counts, "rnaseq-counts.RDS")
write.table(assays(counts)$counts, file="rnaseq-counts.txt")
tsmsg("Counting complete")
