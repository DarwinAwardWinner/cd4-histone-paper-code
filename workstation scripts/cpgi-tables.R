library(Biobase)
library(RColorBrewer)
library(Matrix)
library(rtracklayer)
library(GenomicRanges)
library(magrittr)
library(dplyr)
source("rnaseq-common2.R", chdir=TRUE)

## Load CPGi data and compute overlaps
cpgi <- read.table("CpGi.Ext.bed", sep="\t")
cpgi <- GRanges(seqnames=cpgi[[1]], ranges=IRanges(start=cpgi[[2]]+1, end=cpgi[[3]]), strand="*")

{
    sexp <- readRDS("promoter-group-counts-1kb.RDS")
    mcols(sexp)$CpGi <- overlapsAny(sexp, cpgi)

    desired.ids <- c("UCSCKG", "ENTREZ")

    cpgi.tables <-
        lapply(desired.ids, . %>% {
                   x <- mcols(sexp)[c(., "CpGi")]
                   if (is(x[[1]], "CharacterList"))
                       x[[1]] %<>% { rtracklayer:::pasteCollapse(.) }
                   x %>% as.data.frame %>%
                       setNames(names(x)) %>%
                           group_by_(names(.)[1]) %>%
                               summarize(CpGi=any(CpGi)) %>%
                                   .[order(.[[1]]),] %>%
                                       as.data.frame
               }) %>% setNames(desired.ids)

    write.xlsx.multisheet(cpgi.tables, "results/CpGi_tables_1kb.xlsx", row.names=FALSE)

    ucsc.with.cpgi <- subset(cpgi.tables$UCSCKG, CpGi)$UCSCKG %>% str_split(",") %>% unlist

    add.cpgi.overlap <- function(infile, outfile) {
        infile %>% read.xlsx(sheet=1, colNames=FALSE) %>% .[[1]] %>%
        {
            xs <- str_split(., ",");
            data.frame(PromoterGroup=rep(., elementLengths(xs)),
                       UCSCKG=unlist(xs))
        } %>%
        transform(CpGi = UCSCKG %in% ucsc.with.cpgi) %>%
        group_by(PromoterGroup) %>% summarize(CpGi=any(CpGi)) %>%
        rename(UCSCKG=PromoterGroup) %>% as.data.frame %T>%
        xlsx::write.xlsx(file=outfile, row.names=FALSE)
    }

    add.cpgi.overlap("/home/ryan/Downloads/Naive H3K4me2 alone.xlsx",
                     "/home/ryan/Downloads/Naive H3K4me2 alone CpGi.xlsx")
    add.cpgi.overlap("/home/ryan/Downloads/Memory H3K4me3.xlsx",
                     "/home/ryan/Downloads/Memory H3K4me3 CpGi.xlsx")

}

{
    sexp <- readRDS("promoter-group-counts-2.5kb.RDS")
    mcols(sexp)$CpGi <- overlapsAny(sexp, cpgi)

    desired.ids <- c("UCSCKG", "ENTREZ")
    cpgi.tables <-
        lapply(desired.ids, . %>% {
                   x <- mcols(sexp)[c(., "CpGi")]
                   if (is(x[[1]], "CharacterList"))
                       x[[1]] %<>% { rtracklayer:::pasteCollapse(.) }
                   x %>% as.data.frame %>%
                       setNames(names(x)) %>%
                           group_by_(names(.)[1]) %>%
                               summarize(CpGi=any(CpGi)) %>%
                                   .[order(.[[1]]),] %>%
                                       as.data.frame
               }) %>% setNames(desired.ids)

    write.xlsx.multisheet(cpgi.tables, "results/CpGi_tables_2.5kb.xlsx", row.names=FALSE)
}
