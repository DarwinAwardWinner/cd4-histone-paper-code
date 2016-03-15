#!/usr/bin/env Rscript

source("common.R", chdir=TRUE)
library(GenomicFeatures)
library(dplyr)
library(magrittr)
library(edgeR)
library(DESeq2)
library(openxlsx)
library(ggplot2)
library(stringr)
library(openxlsx)
library(graphite)
library(GSEABase)
library(BiocParallel)
library(foreach)
library(doParallel)

## Can't do parallel because memory :(
## registerDoParallel()
registerDoSEQ()
register(SerialParam())

## Load differential binding data
load("promoter-edger-results3-groups-1kb.RDa")

dge <- edger.analyses$D0.D1$dge
subdesign <- design %>% .[colnames(dge),] %>% .[,colSums(.) > 0]
ct <- "H3K4me3.Naive.D1 - H3K4me3.Naive.D0"
ctmat <- makeContrasts(contrasts=ct, levels=subdesign)

# MSigDB Gene Sets
wanted.gene.sets <-
    c("h", "c1",
      "c2", "c2.CGP",
      "c2.CP", "c2.CP:BIOCARTA", "c2.CP:KEGG", "c2.CP:REACTOME",
      "c3", "c3.MIR", "c3.TFT",
      "c4", "c4.CGN", "c4.CM",
      "c5", "c5.BP", "c5.CC", "c5.MF",
      "c6", "c7")
msigdb <- readRDS("~/Projects/msigdb/msigdb-entrez.RDS") %>%
    .[wanted.gene.sets]
## Convert gene sets to lists of character vectors and subset each
## gene set to only genes represented in the DGEList
msigdb.gsets <- msigdb %>%
    lapply(. %>% geneIds %>% ids2indices(identifiers=dge$genes$ENTREZ)) %>%
    setNames(make.names(names(.)))

## MSigDB CAMERA results
camera.results <- msigdb.gsets %>%
    lapply(camera, y=dge, design=subdesign, contrast=ctmat)
fname <- "results/camera-MSigDB-H3K4me3-N0vN1.xlsx"
write.xlsx.multisheet(camera.results, fname)

## SPIA on graphite pathways

## Prepare the pathways
pathdb.table <- pathwayDatabases() %>%
    lapply(as.character) %>% as_data_frame %>%
    subset(species=="hsapiens") %>%
    mutate(dbname=sprintf("%s_%s", species, database)) %>%
    mutate(SPIAfile=sprintf("%sSPIA.RData", dbname))

if (file.exists("graphite-hsapiens-pathways.RDS")) {
    pathdbs <- readRDS("graphite-hsapiens-pathways.RDS")
} else {
    ## Fetch pathways, convert to Entrez IDs, prepare SPIA files
    pathdbs <- do.call(foreach, pathdb.table) %dopar% {
        db <- pathways(species, database) %>%
            convertIdentifiers("entrez")
        db
    } %>% setNames(pathdb.table$dbname)
    saveRDS(pathdbs, "graphite-hsapiens-pathways.RDS")
}

## These pathways cause prepareSPIA to throw an error
blacklist.pathways <- c(
    "Insulin receptor signalling cascade", "IRS activation",
    "IRS-related events", "SHC activation", "SHC-related events",
    "Signaling by Insulin receptor", "Signaling Pathways"
)
blacklist.regexp <- str_c("(?:", str_c("\\Q_", blacklist.pathways, "\\E", collapse="|"), ")$")

## do.call(foreach, pathdb.table) %dopar% {
##     if (!file.exists(SPIAfile)) {
##         db <- pathdbs[[dbname]] %>%
##             .[! names(.) %in% blacklist.pathways]
##         message("Preparing ", dbname, " for SPIA.")
##         prepareSPIA(db, dbname, TRUE)
##     }
## }

pathdb <- lapply(pathdbs, as.list) %>% unlist
for (i in seq_len(length(pathdb))) {
    pathdb[[i]] %<>% { .@title <- str_c(.@database, "_", .@title); . }
}
names(pathdb) <- vapply(pathdb, . %>% .@title, "")

. %>% { setNames(as.list(.), str_c(.@name, "_", names(.))) }) %>%
    unlist2 %>% new("PathwayList",
                    name="AllDatabases",
                    species=pathdbs[[1]]@species,
                    entries=.,
                    timestamp=pathdbs[[1]]@timestamp)
for(i in names(pathdb)) {
    pathdb[[]]
}
pathdb <- new("PathwayList",
              name="AllDatabases",
              species=pathdbs[[1]]@species,
              entries=pathdb,
              timestamp=pathdbs[[1]]@timestamp)

spiabase <- "graphite"
spiafile <- str_c(spiabase, "SPIA.RData")
if (!file.exists(spiafile))
    prepareSPIA(pathdb %>% .[!str_detect(names(.), blacklist.regexp)],
                spiabase)

## Run SPIA
spia.results <- all.results.ql$H3K4me3.Naive.D0vD1 %$% {
    ENTREZ %<>% as.character
    de <- logFC %>% setNames(ENTREZ) %>% .[FDR <= 0.2] %>%
        ## For ENTREZ IDs with multiple promoters, use the one with
        ## the greatest absolute change.
        .[order(abs(.), decreasing = TRUE)] %>% .[!duplicated(names(.))]
    runSPIA(de=de, all=ENTREZ, spiabase, verbose=TRUE)
}

fname <- "results/SPIA-graphite-H3K4me3-N0vN1.xlsx"
openxlsx::write.xlsx(spia.results, fname)

## Run CAMERA on pathways
pathdb.gset <- pathdb %>% as.list %>% lapply(nodes) %>%
    ids2indices(identifiers=dge$genes$ENTREZ)
graphite.camera.results <- camera(y=dge, index=pathdb.gset, design=subdesign, contrast=ctmat)
fname <- "results/camera-graphite-H3K4me3-N0vN1.xlsx"
openxlsx::write.xlsx(graphite.camera.results, fname)

## Run SPIA on original data
x <- openxlsx::read.xlsx("~/Desktop/H3K4me3 N0vN1.xlsx")
spia.results.ss <- x %$% {
    de <- logFC.H3K4me3 %>% setNames(ENTREZ) %>% .[FDR.H3K4me3 <= 0.1] %>%
        ## For ENTREZ IDs with multiple promoters, use the one with
        ## the greatest absolute change.
        .[order(abs(.), decreasing = TRUE)] %>% .[!duplicated(names(.))]
    ENTREZ <- c(as.character(ENTREZ),
                as.character(all.results.ql$H3K4me3.Naive.D0vD1$ENTREZ)) %>%
        unique
    runSPIA(de=de, all=ENTREZ, spiabase, verbose=TRUE)
}

fname <- "results/SPIA-graphite-H3K4me3-N0vN1-fromSarah.xlsx"
openxlsx::write.xlsx(spia.results.ss, fname)
