library(GenomicRanges)
library(Rsamtools)
library(Rsubread)
library(BatchJobs)
library(BiocParallel)
library(doParallel)
library(parallel)
library(plyr)
library(stringr)
library(gtools)
library(xlsx)
library(BSgenome.Hsapiens.UCSC.hg19)
library(bigmemory)
library(bigmemoryExtras)
library(rtracklayer)

num.cores <- as.numeric(Sys.getenv("PBS_NUM_PPN"))
if (is.na(num.cores))
    num.cores <- parallel:::detectCores()
registerDoParallel(cores=num.cores)
if (num.cores > 1) {
    register(MulticoreParam(num.cores))
} else {
    register(SerialParam())
}
## Not sure which of these is the real option. Some code seems to
## use different ones.
options(cores=num.cores)
options(mc.cores=num.cores)

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

## sprintf with shell-quoting of all substituted arguments
sqprintf <- function(fmt, ...) {
    args <- c(list(fmt=fmt), lapply(list(...), shQuote))
    do.call(sprintf, args)
}

fixdfchar <- function(df) {
    for (i in seq_along(df)) {
        if (is.factor(df[[i]]) && !any(duplicated(df[[i]]))) {
            df[[i]] <- as.character(df[[i]])
        }
    }
    df
}

## Return TRUE if every element of x has a distinct, non-missing name
has.valid.names <- function(x) {
    nm <- names(x)
    !is.null(nm) && !any(is.na(nm)) && !anyDuplicated(nm)
}

## Make seqlengths work on BamFileList
setMethod("seqinfo", signature=list(x="BamFileList"), function (x) {
    Reduce(merge, lapply(x, seqinfo))
})

## Get the chromosome lengths from BSgenome, BamFileList, GRanges, etc.
get.seqlengths <- function(object) {
    sl <- tryCatch(seqlengths(object), error=function(...) {
        setNames(as.numeric(object), names(object))
    })
    stopifnot(has.valid.names(sl))
    sl
}

gr.to.saf <- function(gr) {
    stopifnot(has.valid.names(gr))
    ## We need to go through DataFrame because some elements may be
    ## Rle or other non-primitive vectors
    x <- DataFrame(Chr=seqnames(gr),
                   Start=start(gr),
                   End=end(gr),
                   Strand=strand(gr),
                   GeneID=names(gr))
    as(x, "data.frame")
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

splitByMaxLength <- function(x, maxlength) {
    nchunks <- ceiling(length(x) / maxlength)
    q <- seq(0, 1, length.out=nchunks+1)
    f <- quantcut(1:length(x), q)
    split(x, f)
}

splitRowsByMaxLength <- function(x, maxlength) {
    nchunks <- ceiling(nrow(x) / maxlength)
    q <- seq(0, 1, length.out=nchunks+1)
    f <- quantcut(1:nrow(x), q)
    split.data.frame(x, f)
}

combine.featureCounts.stat <- function(...) {
    x <- list(...)
    stopifnot(length(x) >= 1)
    if (length(x) == 1) return(x)
    values <- lapply(x, function(df) as.matrix(df[-1]))
    combined.values <- Reduce(`+`, values)
    data.frame(x[[1]][1], combined.values)
}

## Wrapper for featureCounts that allows easy specification of
## fragment length
featureCounts.fragments <- function(bam, saf, fraglength, nthreads=max(getOption("cores"), getOption("mc.cores")),
                                    bigmatrix.backing.file, ...) {
    stopifnot(length(fraglength) == 1 && fraglength >= 1)
    fraglength <- floor(fraglength)
    if (is(saf, "GRanges"))
        saf <- gr.to.saf(saf)
    max.chunksize <- min(1e6, floor(.Machine$integer.max / length(bam) / 2))

    safsplit <- unname(splitRowsByMaxLength(saf, max.chunksize))
    if (length(safsplit) > 1) {
        tsmsg("Split long annotation into ", length(safsplit), " chunks of ", max.chunksize ," elements each")
    } else {
        tsmsg("Counting all ", nrow(safsplit[[1]]), " features in one run")
    }
    if (missing(bigmatrix.backing.file)) {
        bigcounts <- big.matrix(
            nrow=nrow(saf), ncol=length(bam), type="integer",
            dimnames=list(GeneID=as.character(saf$GeneID), Sample=names(bam)))
    } else {
        if (file.exists(bigmatrix.backing.file))
            file.remove(bigmatrix.backing.file)
        stopifnot(!file.exists(bigmatrix.backing.file))
        bigcounts <- BigMatrix(
            nrow=nrow(saf), ncol=length(bam), type="integer",
            dimnames=list(GeneID=as.character(saf$GeneID), Sample=names(bam)),
            backingfile=bigmatrix.backing.file)
    }

    resLists <- list(annotation=list(),
                     targets=list(),
                     stat=list())
    i <- 0
    nextrow <- 1
    for (subsaf in safsplit) {
        i <- i+1
        tsmsg("Counting feature chunk ", i)
        x <- featureCounts(
            bam, annot.ext=subsaf, isPairedEnd=FALSE,
            read2pos=5, readExtension3=fraglength-1,
            allowMultiOverlap=TRUE,
            strandSpecific=0, nthreads=nthreads, ...)
        ct <- x$counts
        bigcounts[seq(from=nextrow, by=1, length.out=nrow(ct)),] <- ct
        for (i in names(resLists)) {
            resLists[[i]] <- c(resLists[[i]], unname(x[i]))
        }
        nextrow <- nextrow + nrow(ct)
        tsmsg("Finshed counting feature chunk ", i)
    }
    res <- list(
        counts=bigcounts,
        annotation=do.call(rbind, unname(resLists$annotation)),
        targets=unlist(unname(resLists$targets)),
        stat=do.call(combine.featureCounts.stat, unname(resLists$stat)))
    res
}

featureCounts.fiveprime <- function(bam, saf, ...) {
    featureCounts.fragments(bam, saf, fraglength=1, ...)
}

## Bedtools-based counting functions & support functions
write.bedtools.seqlengths <- function(object, file) {
    object <- get.seqlengths(object)
    object <- data.frame(chrom=names(object), size=object)
    write.table(object, file=file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
}

shared.tempdir <- Sys.getenv("PBSREMOTEDIR")
(function() {
    stopifnot(file.exists(shared.tempdir))
    tf <- tempfile(tmpdir=shared.tempdir)
    writeLines("hello", tf)
    stopifnot(readLines(tf) == "hello")
    TRUE
})()

## Override tempfile to add files to autoclean list
tempfiles.to.delete <- character(0)
tempfile <- function(..., tmpdir=shared.tempdir) {
    ret <- base::tempfile(..., tmpdir=tmpdir)
    tempfiles.to.delete <<- c(tempfiles.to.delete, ret)
    invisible(ret)
}

clean.tempfiles <- function() {
    ret <- logical(0)
    if (length(tempfiles.to.delete))
        suppressWarnings(ret <- file.remove(tempfiles.to.delete))
    tempfiles.to.delete <<- character(0)
    invisible(ret)
}

write.bedtools.seqlengths.to.tempfile <- function(object, ...) {
    file <- tempfile(...)
    write.bedtools.seqlengths(object, file)
    file
}

str_is_whitespace <- function (string) {
    !str_detect(string, perl("\\S"))
}

str_empty <- function(string) {
    is.na(string) | string == ""
}

pipeline <- function(..., infile, outfile) {
    allcmd <- unlist(list(...))
    allcmd <- allcmd[!str_empty(allcmd) & str_detect(allcmd, perl("\\S"))]
    full.pipeline <- str_c(allcmd, collapse=" | ")
    if (!missing(infile))
        full.pipeline <- str_c(sqprintf("<%s ", infile), full.pipeline)
    if (!missing(outfile))
        full.pipeline <- str_c(full.pipeline, sqprintf(" >%s", outfile))
    full.pipeline
}

pipelinevec <- function(...) {
    args <- list(...)
    args <- args[!sapply(args, is.null)]
    args$FUN <- pipeline
    do.call(mapply, args)
}

run.command <- function(cmd, ...) {
    system2("sh", args=c("-c", shQuote(cmd)), ...)
}

## Return pipeline to convert ranges to just their 5-prime ends (must
## be stranded, result is zero-length ranges). Minor bug: will never
## return a 5-prime end at the first or last base of a chromosome.
bedtools.length.to.zero.cmd <- function() {
    ## This perl script does the same as the below
    "perl -lane 'if ($F[5] eq q(+)) { $F[2] = $F[1] } else { $F[1] = $F[2] }; print (join(q(\t), @F))"
    ## pipeline(
    ##     sqprintf(c("bedtools flank -i /dev/stdin -l 1 -r 0 -s -g %s",
    ##                "bedtools flank -i /dev/stdin -l 0 -r 1 -s -g %s",
    ##                "bedtools slop -i /dev/stdin -l 0 -r -1 -s -g %s"),
    ##              genome.file))
}

## Adjust range length to specified length, anchoring 5-prime end and
## ignoring original 3-prime end.
bedtools.setlength.cmd <- function(length) {
    length <- as.numeric(length)
    ## This perl script does the same as the below
    sqprintf("perl -MList::Util -lane 'if ($F[5] eq q(+)) { $F[2] = $F[1] + %s } else { $F[1] = $F[2] - %s }; $F[1] = List::Util::max($F[1], 0); print (join(q(\t), @F))'", length, length)
    ## pipeline(
    ##     bedtools.length.to.zero.cmd(genome.file),
    ##     sqprintf("bedtools slop -i /dev/stdin -l 0 -r %s -s -g %s",
    ##              length, genome.file))
}

bedtools.shift.cmd <- function(shift) {
    shift <- as.numeric(shift)
    if (shift != 0) {
        ## This perl script does the same as the below
        sqprintf("perl -lane 'if ($F[5] eq q(+)) { $F[1] = $F[1] + %s; $F[2] = $F[2] + %s } else { $F[1] = $F[1] - %s; $F[2] = $F[2] - %s }; $F[1] = List::Util::max($F[1], 0); $F[2] = List::Util::max($F[2], 0); print (join(q(\t), @F))'", length, length, length, length)

        ## sqprintf("bedtools slop -i /dev/stdin -l %s -r %s -s -g %s",
        ##          -length, length, genome.file)
    }
}

countFragmentOverlapsByBedtools <- function(bam, annot, fraglength, shift=0, bigmatrix.backing.file, ...) {
    tsmsg("Checking inputs")
    on.exit(clean.tempfiles())
    stopifnot(length(fraglength) == 1 && fraglength >= 1)
    fraglength <- floor(fraglength)
    stopifnot(has.valid.names(annot))
    tsmsg("Setting up bigmatrix")
    if (missing(bigmatrix.backing.file)) {
        bigcounts <- big.matrix(
            nrow=length(annot), ncol=length(bam), type="integer",
            dimnames=list(Feature=as.character(names(annot)), Sample=names(bam)))
    } else {
        if (file.exists(bigmatrix.backing.file))
            file.remove(bigmatrix.backing.file)
        stopifnot(!file.exists(bigmatrix.backing.file))
        bigcounts <- BigMatrix(
            nrow=length(annot), ncol=length(bam), type="integer",
            dimnames=list(Feature=as.character(names(annot)), Sample=names(bam)),
            backingfile=bigmatrix.backing.file)
    }
    tsmsg("Writing annotation to BED file")
    annotfile <- tempfile("annot", fileext=".bed")
    annot$name <- names(annot)
    export(annot, annotfile, format="bed")

    tsmsg("Constructing bedtools counting pipelines")
    output.tempfiles <- tempfile(str_c("output", seq_along(bam)))
    cmds <- pipelinevec(
        sqprintf("bedtools bamtobed -i %s", bam),
        ## Set fraglength
        bedtools.setlength.cmd(fraglength),
        ## Apply shift
        bedtools.shift.cmd(shift),
        ## Compute coverage counts
        sqprintf("bedtools coverage -a /dev/stdin -b %s -counts", annotfile),
        ## Extract only feature name & count
        sqprintf("awk '{print $4,$(NF)}'"),
        outfile=output.tempfiles)

    tsmsg("Executing bedtools counting pipelines")
    bp <- BatchJobsParam(resources=list(modules=c("bedtools"), nodes="1:ppn=2"))
    bplapply(cmds, system, BPPARAM=bp)
    tsmsg("Loading counts into bigmatrix")
    bplapply(seq_along(bam), function(i) {
        result.table <- read.table(output.tempfiles[i], as.is=TRUE, header=F)
        names(result.table) <- c("Feature", "Count")
        bigcounts[result.table$Feature, i] <- result.table$Count
    })
    tsmsg("Finished counting")
    bigcounts
}

data.frame2GRanges <- function(df, keepColumns = FALSE, ignoreStrand = FALSE,
                               startOffset=0, endOffset=0) {
    stopifnot(class(df) == "data.frame")
    stopifnot(all(c("start", "end") %in% names(df)))
    stopifnot(any(c("chr", "seqnames") %in% names(df)))
    if("seqnames" %in% names(df))
        names(df)[names(df) == "seqnames"] <- "chr"
    if(!ignoreStrand && "strand" %in% names(df)) {
        if(is.numeric(df$strand)) {
            strand <- ifelse(df$strand == 1, "+", "*")
            strand[df$strand == -1] <- "-"
            df$strand <- strand
        }
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start + startOffset, end = df$end + endOffset),
                      strand = df$strand)
    } else {
        gr <- GRanges(seqnames = df$chr,
                      ranges = IRanges(start = df$start + startOffset, end = df$end + endOffset))
    }
    if(keepColumns) {
        dt <- as(df[, setdiff(names(df), c("chr", "start", "end", "strand"))],
                 "DataFrame")
        elementMetadata(gr) <- dt
    }
    names(gr) <- rownames(df)
    gr
}

## Return base positions separated by the specified interval across a
## chromosome such that the positions are symmetric on both ends (plus
## or minus 1 bp if the length is odd). These are suitable for use as
## the centers of sliding windows on that chromosome.
symmetric.intervals <- function(seqlength, interval) {
    seq(from=1, to=seqlength, by=interval) + floor(((seqlength - 1) %% interval) / 2)
}

## Generate sliding windows of specified size at specified intervals
## along each chromsome in genome. Genome can be a n
sliding.windows <- function(genome, window.interval, window.size=window.interval) {
    sl <- get.seqlengths(genome)
    windows <- IRangesList(bplapply(sl, function(l) {
        i <- symmetric.intervals(l, window.interval)
        ir <- IRanges(i - floor( (window.size-1) / 2), i + ceiling( (window.size-1) / 2))
        start(ir) <- pmax(start(ir), 1)
        end(ir) <- pmin(end(ir), l)
        ir
    }))
    x <- GRanges(seqnames=Rle(names(windows), elementLengths(windows)),
                 ranges=unlist(unname(windows)), strand=Rle("*"),
                 seqlengths=sl)
    names(x) <- sprintf("%s:%s-%s", as.character(seqnames(x)), start(x), end(x))
    x
}

read.narrowPeak <- function(file, ...) {
    peaks.df <- read.table(file, sep="\t", row.names=NULL, ...)
    names(peaks.df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
    peaks.df$strand <- "*"
    peaks.df$name <- as.character(peaks.df$name)
    row.names(peaks.df) <- str_replace(as.character(peaks.df$name), "^.*_peak", "peak")
    data.frame2GRanges(peaks.df, keepColumns=TRUE, startOffset=1, endOffset=0)
}

write.narrowPeak <- function(x, file, ...) {
    x <- as(x, "data.frame")
    if("seqnames" %in% names(x))
        names(x)[names(x) == "seqnames"] <- "chr"
    x <- x[c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")]
    write.table(x, file, sep="\t", row.names=FALSE, col.names=FALSE, ...)
}
