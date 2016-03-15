#!/usr/bin/Rscript

source("rnaseq-common2.R", chdir=TRUE)

library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(annotate)
library(rtracklayer)

flatten.CharacterList <- function(x, sep=",") {
    tryCatch(x <- as(x, "CharacterList"),
             error=function(...) NULL)
    if (is(x, "CharacterList"))
        setNames(rtracklayer:::pasteCollapse(x, collapse=sep),
                 names(x))
    else
        x
}

fix.cl.in.df <- function(df) {
    df[] <- llply(df, flatten.CharacterList, .parallel=TRUE)
}

load.object.from.file <- function(datafile, varname) {
    within(list(), load(datafile))[[varname]]
}

alltestnames <- c("Naive.AllT", "Memory.AllT", "Naive.D0vD1",
                  "Naive.D0vD5", "Naive.D0vD14", "Memory.D0vD1",
                  "Memory.D0vD5", "Memory.D0vD14", "NvM.AllT", "NvM.D0",
                  "NvM.D1", "NvM.D5", "NvM.D14", "Fac.AllT", "Fac.D1",
                  "Fac.D5", "Fac.D14", "MD0vND14")

{
    tsmsg("Loading 1kb promoter data")
    promoter.1kb.alltabs <- local({
        edger.analyses <- load.object.from.file("promoter-edger-results3-groups-1kb.RDa", "edger.analyses")
        tsmsg("Extracting promoter result tables")
        all.results.ql <- do.call(c, unname(plyr::llply(edger.analyses, `[[`, "results.ql")))
        rm(edger.analyses)
        gc()
        all.results.ql <- plyr::llply(all.results.ql, function(x) {
            x$Input.signif <- rownames(x) %in% rownames(all.results.ql$input.consistency)
            reorder.columns(x, start=c("ENTREZ", "SYMBOL", "GENENAME"))
        })
        names(all.results.ql) <- str_replace(names(all.results.ql), "M0vN14", "MD0vND14")
        all.results.ql[str_detect(names(all.results.ql), "H3K4me[23]")]
    })
    tsmsg("Loading 2.5kb promoter data")
    promoter.2.5kb.alltabs <- load.object.from.file("promoter-edger-results3-groups-2.5kb.RDa", "edger.analysis")$results.ql
    tsmsg("Loading RNA-seq tables")
    rnaseq.alltabs <- load.object.from.file("rnaseq-edger-results3.RDa", "edger.analysis")$results.ql

    tsmsg("Grouping tables by sample type")
    me2tabs <- llply(alltestnames, function(i) promoter.1kb.alltabs[[str_c("H3K4me2.", i)]])
    me3tabs <- llply(alltestnames, function(i) promoter.1kb.alltabs[[str_c("H3K4me3.", i)]])
    k27tabs <- llply(alltestnames, function(i) promoter.2.5kb.alltabs[[str_c("H3K27me3.", i)]])
    rnatabs <- rnaseq.alltabs[alltestnames]
    stopifnot(!any(sapply(list(me2tabs, me3tabs, k27tabs, rnatabs), is.null)))
    stopifnot(!any(sapply(c(me2tabs, me3tabs, k27tabs, rnatabs), is.null)))

    rm(promoter.1kb.alltabs, promoter.2.5kb.alltabs, rnaseq.alltabs); gc()
}

pre.metacols <- c("ENTREZ", "SYMBOL", "GENENAME", "UCSCKG", "CpGi")
post.metacols <- c("UNIGENE", "ENSEMBL", "REFSEQ", "UNIPROT")
statcols <- c("PValue", "FDR", "logFC")

## Makre sure same TSS list
allkg.1kb <- unlist(str_split(me2tabs[[1]]$UCSCKG, fixed(",")))
allkg.2.5kb <- unlist(str_split(k27tabs[[1]]$UCSCKG, fixed(",")))
stopifnot(length(allkg.2.5kb) == length(allkg.1kb))
stopifnot(length(intersect(allkg.1kb, allkg.2.5kb)) == length(allkg.2.5kb))
stopifnot(length(setdiff(allkg.1kb, allkg.2.5kb)) + length(setdiff(allkg.2.5kb, allkg.1kb)) == 0)

mergecol.2.5kb <- seq(nrow(k27tabs[[1]]))
mergecol.1kb <- local({
    first.kg.1kb <- str_extract(me2tabs[[1]]$UCSCKG, "^[^,]*")
    all.kg.2.5kb <- k27tabs[[1]]$UCSCKG
    laply(llply(first.kg.1kb, grep, all.kg.2.5kb, fixed=TRUE), `[`, 1)
})

tsmsg("Extracting stat columns")
promoter.1kb.stat.tabs <- llply(alltestnames, function(i) {
    x <- list(H3K4me2=me2tabs[[i]],
              H3K4me3=me3tabs[[i]])
    x <- llply(x, function(x) x[na.omit(match(statcols, colnames(x)))])
    x <- llply(names(x), function(nm) setNames(x[[nm]], str_c(names(x[[nm]]), ".", nm)))
    do.call(data.frame, unname(x))
})
promoter.2.5kb.stat.tabs <- llply(alltestnames, function(i) {
    x <- list(H3K27me3=k27tabs[[i]])
    x <- llply(x, function(x) x[na.omit(match(statcols, colnames(x)))])
    x <- llply(names(x), function(nm) setNames(x[[nm]], str_c(names(x[[nm]]), ".", nm)))
    do.call(data.frame, unname(x))
})
promoter.stat.tabs <- llply(names(promoter.1kb.stat.tabs), function(i) {
    df.1kb <- data.frame(mergecol=mergecol.1kb, promoter.1kb.stat.tabs[[i]])
    df.2.5kb <- data.frame(mergecol=mergecol.2.5kb, promoter.2.5kb.stat.tabs[[i]])
    merge(df.1kb, df.2.5kb, by="mergecol", all=TRUE)[-1]
})
names(promoter.stat.tabs) <- names(promoter.1kb.stat.tabs)
rna.stat.tabs <- llply(rnatabs, function(x) {
    x <- x[na.omit(match(statcols, colnames(x)))]
    names(x) <- str_c(names(x), ".RNASeq")
    x
})

## pre.meta.tab <- gene.annot[pre.metacols]
## post.meta.tab <- gene.annot[post.metacols]

pre.meta.tab <- me2tabs[[1]][pre.metacols]
post.meta.tab <- me2tabs[[1]][post.metacols]

## Find matching RNA result for each
rna.matched.rows <- match(as.character(me2tabs[[1]]$ENTREZ), as.character(rnatabs[[1]]$ENTREZ))

## all.rownames <- sort(union(rownames(promoter.stat.tabs[[1]]), rownames(rna.stat.tabs[[1]])))

tsmsg("Merging all tables")
full.tabs <- llply(alltestnames, function(i) {
    data.frame(pre.meta.tab,
               promoter.stat.tabs[[i]],
               rna.stat.tabs[[i]][rna.matched.rows,],
               post.meta.tab)
}, .parallel=TRUE)

tsmsg("Saving merged tables")
saveRDS(full.tabs, "edger-merged-tables3-ql.RDS")

toptabs <- llply(full.tabs, function(tab) {
    bestfdr <- do.call(pmin, c(na.rm=TRUE, tab[str_detect(names(tab), "^FDR")]))
    tab <- tab[bestfdr <= 0.1,]
    bestfdr <- bestfdr[bestfdr <= 0.1]
    tab <- tab[order(bestfdr),]
    tab
}, .parallel=TRUE)

write.xlsx.multisheet(toptabs, "results/edger-merged-tables3-ql.xlsx", row.names=FALSE)
