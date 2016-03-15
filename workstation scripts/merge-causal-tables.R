#!/usr/bin/Rscript

source("rnaseq-common.R", chdir=TRUE)

promoter <- within(list(), load("promoter-edger.RDa"))[c("timepoint.analysis", "celltype.analysis", "factorial.analysis")]
rnaseq <- within(list(), load("rnaseq-edger.RDa"))[c("timepoint.analysis", "celltype.analysis", "factorial.analysis")]

nametrans <- c(Naive="N", Memory="M")

out.tables <- list()

## Get ChIP N0 vs N120 and RNA-seq N0 vs M0 in one table
{
    di <- promoter$timepoint.analysis$IP.diH3K4.N$results.lr
    tri <- promoter$timepoint.analysis$IP.triH3K4.N$results.lr
    rna <- rnaseq$celltype.analysis$T0$results

    ## Pick best TSS PValue for each gene
    di <- di[order(di$PValue),]
    tri <- tri[order(tri$PValue),]
    di <- di[!duplicated(di$geneSymbol),]
    tri <- tri[!duplicated(tri$geneSymbol),]

    rna <- rna[rna$SYMBOL != "",]

    di <- with(di,
               data.frame(row.names=kgID,
                          SYMBOL=geneSymbol,
                          diH3K4.PValue=PValue,
                          diH3K4.FDR=FDR,
                          diH3K4.logFC=logFC))
    tri <- with(tri,
                data.frame(row.names=kgID,
                           SYMBOL=geneSymbol,
                           triH3K4.PValue=PValue,
                           triH3K4.FDR=FDR,
                           triH3K4.logFC=logFC))
    allprom <- merge(di, tri, by="SYMBOL", sort=FALSE)
    rna <- rna[c("PValue", "FDR", "logFC", "ENTREZ", "SYMBOL", "GENENAME", "UNIGENE", "ENSEMBL", "REFSEQ", "UCSCKG", "UNIPROT")]
    rna <- rename(rna, c(PValue="rna.PValue",
                         FDR="rna.FDR",
                         logFC="rna.logFC"))
    alldata <- merge(allprom, rna, by="SYMBOL", all=TRUE)
    out.tables[["ChIP T0vT120 & RNA NvM"]] <- alldata
}

## Get ChIP N0 vs M0 and RNA-seq N0 vs N120 in one table
{
    di <- promoter$celltype.analysis$IP.diH3K4.T0$results.lr
    tri <- promoter$celltype.analysis$IP.triH3K4.T0$results.lr
    rna <- rnaseq$timepoint.analysis$Naive$results

    ## Pick best TSS PValue for each gene
    di <- di[order(di$PValue),]
    tri <- tri[order(tri$PValue),]
    di <- di[!duplicated(di$geneSymbol),]
    tri <- tri[!duplicated(tri$geneSymbol),]

    rna <- rna[rna$SYMBOL != "",]

    di <- with(di,
               data.frame(row.names=kgID,
                          SYMBOL=geneSymbol,
                          diH3K4.PValue=PValue,
                          diH3K4.FDR=FDR,
                          diH3K4.logFC=logFC))
    tri <- with(tri,
                data.frame(row.names=kgID,
                           SYMBOL=geneSymbol,
                           triH3K4.PValue=PValue,
                           triH3K4.FDR=FDR,
                           triH3K4.logFC=logFC))
    allprom <- merge(di, tri, by="SYMBOL", sort=FALSE)
    rna <- rna[c("PValue", "FDR", "logFC", "ENTREZ", "SYMBOL", "GENENAME", "UNIGENE", "ENSEMBL", "REFSEQ", "UCSCKG", "UNIPROT")]
    rna <- rename(rna, c(PValue="rna.PValue",
                         FDR="rna.FDR",
                         logFC="rna.logFC"))
    alldata <- merge(allprom, rna, by="SYMBOL", all=TRUE)
    out.tables[["ChIP NvM & RNA T0vT120"]] <- alldata
}

out.tables.signif <- llply(out.tables, function(x) {
    signif <- which(with(x, pmin(diH3K4.FDR, triH3K4.FDR, rna.FDR) <= 0.1))
    x <- x[signif,]
    x[order(with(x, pmin(diH3K4.FDR, triH3K4.FDR, rna.FDR))),]
})

write.xlsx.multisheet(out.tables.signif, "edger-causal-tables.xlsx")
