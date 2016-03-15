suppressMessages(suppressWarnings({
  ## Prevent gsubfn from requiring tcl (at the cost of a speed decrease)
  options(gsubfn.engine = "R")

  ## Make sure we start Java with lots of heap space before anything
  ## else triggers it to start with the defaults
  options(java.parameters=c("-Xrs", "-Xmx8g"))
  library(stringr)
  library(foreach)
  library(parallel)
  library(doParallel)
  ## Not sure which of these is the real option. Some code seems to
  ## use different ones.
  options(cores=parallel:::detectCores())
  options(mc.cores=parallel:::detectCores())
  options(mc.preschedule=FALSE)
  registerDoParallel()
  library(plyr)
  library(Rsamtools)
  library(ShortRead)
  library(GenomicRanges)
  library(Biobase)
  library(Biostrings)
  library(rtracklayer)
  library(affy)
  library(ggplot2)
  library(inline)
  library(stringr)
  library(chipseq)
  library(GenomicFeatures)
  ## library(spp)
  library(plyr)
  ## library(R.utils)
  library(sqldf)
  ## library(openxlsx)
  library(functional)
  library(biomaRt)
}))

source("fixes.R")

renice.self <- function() {
  pid <- Sys.getpid()
  system2("ionice", c("-c3", "-p", pid), stdout=FALSE)
  system2("renice", c("-n", 19, "-p", pid), stdout=FALSE)
}

makeNiceCluster <- function(...) {
  cl <- makeCluster(...)
  clusterCall(cl, renice.self)
  cl
}

tsmsg <- function(...) {
  message(date(), ": ", ...)
}

tsmsgf <- function(...) {
  tsmsg(sprintf(...))
}

multiplicity <- function(v) {
  table(v)[as(v, "vector")]
}

seqlengthsForUCSCGenome <- function(genome) {
  seqlengths(GRangesForUCSCGenome(genome))
}

chromRangesLists <- function(chromosomes) {
  chr.ranges <- IRanges(1, chromosomes, names=names(chromosomes))
  chr.rangeslists <- foreach(i=1:length(chromosomes)) %do% {
    rl <- RangesList(chr.ranges[i]);
    names(rl) <- names(chromosomes)[i]
    rl
  }
  names(chr.rangeslists) <- names(chromosomes)
  return(chr.rangeslists)
}

readAlignedRanges <- function(bamfile, genome=NULL, include.reads=TRUE) {
  bamheader <- scanBamHeader(bamfile)[[1]]
  message("Reading ", bamfile)
  ## Generate a separate RangesList for each chromosome, so that we
  ## can load reads from one chromosome at a time to save memory.
  chr.rangeslists <- chromRangesLists(bamheader$targets)
  ## Load the reads from the BAM file one chromosome at a time
  full.gr <- foreach(chr.rl=chr.rangeslists, .combine=c) %dopar% {
    message("Reading chromosome ", names(chr.rl)[[1]])
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE), which=chr.rl,
                          simpleCigar = TRUE, reverseComplement = TRUE, what = ShortRead:::.readAligned_bamWhat())
    reads <- readAligned(bamfile, type="BAM", param=param)
    gr <- as(reads, "GRanges")
    if (include.reads) {
      elementMetadata(gr)$seq <- sread(reads)
      elementMetadata(gr)$qual <- as(quality(reads), "PhredQuality")
    }
    rm(reads)
    ## Convert flag to Rle so it takes up less memory
    elementMetadata(gr)$flag <- Rle(elementMetadata(gr)$flag)
    message("Finished reading chromosome ", names(chr.rl)[[1]])
    gr
  }
  if (!is.null(genome)) {
    ## If genome is given, just verify that the seqlengths match those given in the bam file.

    x <- seqlengths(GRangesForUCSCGenome(genome))
    y <- bamheader$targets
    stopifnot(all(names(y) %in% names(x)))
    stopifnot(all(y == x[names(y)]))
  } else {
    ## No genome provide, so just use the seqlengths from the BAM file
    seqlengths(full.gr) <- bamheader$targets
  }
  message("Finished reading ", bamfile)
  gc()
  full.gr
}

## Eval expression in forked process in the foreground (not in parallel)
in.forked.process <- function(expr) {
  result <- mccollect(mcparallel(expr))
  if (length(result) == 0) {
    NULL
  } else {
    result <- result[[1]]
    if (is(result, "try-error")) {
      stop(result)
    } else {
      result
    }
  }
}

## Quarantine xlsx package to subprocess
read.xlsx <- function(...) {
    args <- list(...)
    in.forked.process({
        library(xlsx)
        do.call(xlsx::read.xlsx, args)
    })
}

## Quarantine xlsx package to subprocess
write.xlsx <- function(...) {
    args <- list(...)
    in.forked.process({
        library(xlsx)
        do.call(xlsx::write.xlsx, args)
    })
}

## Additional args are passed to every call of addDataFrame
write.xlsx.multisheet <- function(data.frames, file, sheetNames=names(data.frames), ...) {
    moreArgs <- list(...)
    in.forked.process({
        library(xlsx)
        if (is.null(sheetNames)) {
            sheetNames <- str_c("Sheet", seq_along(data.frames))
        }
        ## Ensure correct number of sheetNames
        stopifnot(length(sheetNames) == length(data.frames))
        ## Fill in missing names if needed
        sheetNames[is.na(sheetNames)] <- str_c("Sheet", seq_along(data.frames))[is.na(sheetNames)]
        wb <- createWorkbook()
        sheets <- llply(sheetNames, function(x) createSheet(wb, sheetName=x))
        plyargs <- c(list(cbind(sheet=sheets, x=data.frames), .fun=addDataFrame, .parallel=FALSE), moreArgs)
        do.call(plyr::mlply, plyargs)
        saveWorkbook(wb, file)
    })
}

## Reading via sqldf is much faster for very large tables
read.table.via.sqldf <- function(filename, header=FALSE, row.names=FALSE, sep="\t", ...) {
  f <- file(filename)
  df <- sqldf("select * from f", dbname=tempfile(), file.format=list(header = header, row.names = row.names, sep=sep))
  close(f)
  df
}

readBedRanges <- function(bedfile, genome=NULL, include.reads=FALSE) {
  df <- read.table.via.sqldf(bedfile)
  x <- GRanges(seqnames = Rle(factor(df[[1]])),
               ranges =
               IRanges(start=as.integer(df[[2]]),
                       end=as.integer(df[[3]])),
               id = BStringSet(as.character(df[[4]])),
               score = Rle(as.integer(df[[5]])),
               strand = Rle(factor(df[[6]])))
  ## Pick up seqlengths from given genome if specified
  if (!is.null(genome)) {
    seqlengths(x) <- seqlengths(GRangesForUCSCGenome(genome))[names(seqlengths(x))]
  }
  rm(df)
  gc()
  x
}

strandedCoverageAsGRanges <- function(gr) {
  gr.byStrand <- split(gr, strand(gr), drop=TRUE)
  cov.byStrand <- Map(function(x) as(coverage(x), "GRanges"), gr.byStrand)
  ## The coverage function is unstranded, so we have to manually set the strand in each GRanges object after computing the coverage
  foreach(s=names(cov.byStrand)) %do% {
    strand(cov.byStrand[[s]]) <- Rle(s)
    seqlengths(cov.byStrand[[s]]) <- seqlengths(gr)
  }
  Reduce(c, cov.byStrand)
}

## Returns a DataFrame with columns "ref", "pos", and "strand",
## representing the reference sequence, start position on that
## reference, and strand ("+" or "-") on that reference where the
## 5-prime end of the read starts. For unstranded rows, pos will be
## NA.
getStrandAndPos <- function(gr) {
  strandData <- strand(gr)
  posData <- ifelse(strandData == "+", start(gr), end(gr))
  unstranded <- strandData == "*"
  if (any(unstranded)) {
    posData[strandData == "*"] <- NA
  }
  DataFrame(ref=seqnames(gr), pos=posData, strand=strandData)
}

.pairwiseSelfDistance.internal.primitive <- {
  sig <- signature(pos_count="integer",
                   pos_vector="integer",
                   limit="integer",
                   result="integer")
  code <- '
for ( int i = 0; i < (*pos_count) - 1; i++ ) {
    for( int j = i + 1; j < *pos_count; j++ ) {
        int d = pos_vector[j] - pos_vector[i];
        if( d > *limit ) break;
        if( d > 0 ) result[d-1]++;
    }
}
'
  cfunction( sig, code, language="C", convention=".C" )
}

.pairwiseSelfDistance.internal.helper <- function(pos.vector, limit) {
  pos.vector <- sort(as.integer(as.vector(pos.vector)))
  result <- rep(0,limit)
  .pairwiseSelfDistance.internal.primitive(length(pos.vector), pos.vector, limit, result)$result
}

pairwiseSelfDistanceTable <- function(sap, limit=500) {
  sap.byref <- split(sap[c("pos", "strand")], sap$ref)
  full.values <- foreach(sap=sap.byref, .combine=`+`) %dopar% {
    pos.bystrand <- split(sap$pos, sap$strand)
    plus.values <- .pairwiseSelfDistance.internal.helper(pos.bystrand$`+`, limit)
    minus.values <- .pairwiseSelfDistance.internal.helper(pos.bystrand$`-`, limit)
    values <- c(rev(minus.values),
                mean(c(plus.values[1],minus.values[1])),
                plus.values)
    values
  }
  data.frame(distance=seq(-limit,limit), count=full.values)
}

.pairwiseInterDistance.internal.primitive <- {
  sig <- signature(plus_count="integer",
                   plus_vector="integer",
                   minus_count="integer",
                   minus_vector="integer",
                   limit="integer",
                   result="integer")
  code <- '
int mstart = 0; /* Minus start position */
for ( int ppos = 0; ppos < *plus_count; ppos++ ) {
    /* Advance mstart to first in-range element of minus_vector */
    while(mstart < *minus_count && plus_vector[ppos] - minus_vector[mstart] > *limit) {
        mstart++;
    }
    if (mstart == *minus_count) {
        /* Reached end of minus_vector */
        break;
    }
    /* Compute differences until we hit the first out-of-range
     * element of minus_vector */
    for ( int mpos = mstart; mpos < *minus_count; mpos++ ) {
        int d = minus_vector[mpos] - plus_vector[ppos];
        if ( d > *limit ) break;
        /* limit is the offset for indexing into result, since result
         * ranges from -limit to +limit. */
        result[d+*limit]++;
    }
}'
  cfunction( sig, code, language="C", convention=".C" )
}

.pairwiseInterDistance.internal.helper <- function(plus.pos.vector, minus.pos.vector, limit) {
  plus.pos.vector <- sort(as.integer(as.vector(plus.pos.vector)))
  minus.pos.vector <- sort(as.integer(as.vector(minus.pos.vector)))
  result <- rep(0,limit * 2 + 1)
  .pairwiseInterDistance.internal.primitive(length(plus.pos.vector), plus.pos.vector,
                                            length(minus.pos.vector), minus.pos.vector,
                                            limit, result)$result
}

pairwiseInterDistanceTable <- function(sap, limit=500) {
  sap.byref <- split(sap[c("pos", "strand")], sap$ref)
  full.values <- foreach(sap=sap.byref, .combine=`+`) %dopar% {
    pos.bystrand <- split(sap$pos, sap$strand)
    values <- .pairwiseInterDistance.internal.helper(pos.bystrand$`+`, pos.bystrand$`-`, limit)
    values
  }
  data.frame(distance=seq(-limit,limit), count=full.values)
}

computeTagAutocorrelation <- function(sap, ...) {
  same.strand.df <- data.frame(pairwiseSelfDistanceTable(sap, ...),
                               strand=factor("same", levels=c("same", "opposite")))
  opposite.strand.df <- data.frame(pairwiseInterDistanceTable(sap, ...),
                                   strand=factor("opposite", levels=c("same", "opposite")))
  rbind(same.strand.df, opposite.strand.df)
}

computeClonalityTable <- function(sap) {
  x <- as.data.frame(table(runLength(as(do.call(paste, args=c(as.list(sap), list(sep="_"))), "Rle")),
                           dnn="Count"))
  x$Count <- as.numeric(as.character(x$Count))
  x
}

file.extension <- function(filename) {
  str_extract(filename, "[^.]*?$")
}

file.readers <- list(bam=readAlignedRanges,
                     bed=readBedRanges)

select.reader <- function (filename) {
  reader <- file.readers[[file.extension(filename)]]
  stopifnot(is.function(reader))
  reader
}

load.chipseq.file <- function(f, genome=NULL) {
  gc()
  message("Working on ", f)
  read.ranges <- select.reader(f)(f, genome=genome, include.reads=FALSE)
  message("Computing autocorrelation & clonality")
  strand.and.pos <- getStrandAndPos(read.ranges)
  autocorr <- computeTagAutocorrelation(strand.and.pos)
  clonality <- computeClonalityTable(strand.and.pos)
  message("Finished with ", f)
  list(reads=read.ranges,
       autocorrelation=autocorr,
       clonality=clonality)
}

load.files <- function(files, genome="hg19", parallel=FALSE) {
  `%do.maybe.par%` <- if (parallel) `%dopar%` else `%do%`
  results <- Map(function(f) load.chipseq.file(f, genome), files)
  names(results) <- basename(infiles)
  results
}

plot.autocorr <- function(ac) {
  ggplot(ac, aes(x=distance, y=count, color=strand)) +
    geom_line() + facet_wrap(~ sample)
}

## Plotting clonality
plot.clonality <- function(clon) {
  ggplot(clon, aes(x=Count, y=Freq)) +
    ## Can choose bar graph or line graph
    geom_bar(aes(x=factor(Count))) + scale_x_discrete("Tag clonality") +
      ## geom_line() + geom_point() + scale_x_continuous("Tag clonality") +
      scale_y_continuous("Count of unique tags") +
        facet_wrap(~ sample)
}


## Functions for retrieving a track or table from UCSC, caching it
## locally
(function() {
  ## Inside this function is a "closure", which effectively gives a
  ## place for private variables, which are used for caching without
  ## polluting the global environment.

  session <- NULL
  ## Wrapper function that autospawns the session the first time it is
  ## asked for.
  setup.session <- function() {
    if (is.null(session)) {
      session <<- browserSession("UCSC")
    }
    stopifnot(!is.null(session))
    return(session)
  }

  set.ucsc.genome <- function(genome=NULL) {
    setup.session()
    if (!is.null(genome)) {
      rtracklayer::genome(session) <- genome
    }
    rtracklayer::genome(session)
  }

  get.ucsc.track <- function (trackname, genome=NULL, cache.dir="ucsc_cache") {
    setup.session()
    genome <- set.ucsc.genome(genome)
    track.filename <- file.path(cache.dir, sprintf("%s.%s.bed", genome, trackname))
    if (file.exists(track.filename)) {
      ## Load cached track file
      track.table <- import(track.filename, format="bed")
    } else {
      ## Download new track file from UCSC
      track.table <- track(session, trackname, range=genome)
      ## Cache the track
      dir.create(path=cache.dir, recursive=TRUE, showWarnings=FALSE)
      export(track.table, track.filename, format="bed")
    }
    return(track.table)
  }

  get.ucsc.table <- function(trackname, tablename, genome=NULL, cache.dir="ucsc_cache") {
    setup.session()
    genome <- set.ucsc.genome(genome)
    table.filename <- file.path(cache.dir, sprintf("%s.%s.%s.tsv", genome, trackname, tablename))
    query <- ucscTableQuery(session, trackname, table=tablename)
    if (file.exists(table.filename)) {
      ## Get the table schema so we can provide a colClasses argument
      ## to read.table
      schema <- ucscSchema(query)
      colClasses <- as.character(schema$RType)
      ## Load cached table file
      table <- read.table(table.filename, sep="\t",
                          header=TRUE, row.names=NULL,
                          colClasses=colClasses)
    } else {
      ## Download new track file from UCSC
      table <- getTable(query)
      ## Replace empty strings in factors and character with NA
      table <- data.frame(lapply(table, function (column) {
        if (is.factor(column)) {
          levels(column)[levels(column) == ""] <- NA
        } else if (is.character(column)) {
          column[column == ""] <- NA
        }
        column
      }))

      ## Cache the table
      dir.create(path=cache.dir, recursive=TRUE, showWarnings=FALSE)
      write.table(table, file=table.filename, sep="\t",
                  row.names=FALSE, col.names=TRUE)
    }
    return(table)
  }

  ## Inject public functions into parent namespace
  foreach(x=c("set.ucsc.genome", "get.ucsc.track", "get.ucsc.table")) %do% {
    assign(x, get(x), pos=sys.frame(0))
  }
})()

## Get a transcript db from UCSC (GenomicFeatures package), caching it locally
get.ucsc.transcript.db <- function (genome="hg19", tablename="knownGene", cache.dir="ucsc_cache") {
  filename <- file.path(cache.dir, sprintf("UCSC.%s.%s.tdb.sqlite", genome, tablename))
  if (file.exists(filename)) {
    tdb <- loadDb(filename)
  } else {
    tdb <- makeTranscriptDbFromUCSC(genome=genome, tablename=tablename)
    dir.create(path=cache.dir, recursive=TRUE, showWarnings=FALSE)
    saveDb(tdb, file=filename)
  }
  tdb
}

get.biomart.transcript.db <- function (biomart="ensembl", dataset="hsapiens_gene_ensembl", cache.dir="ucsc_cache") {
  filename <- file.path(cache.dir, sprintf("Biomart.%s.%s.tdb.sqlite", biomart, dataset))
  if (file.exists(filename)) {
    tdb <- loadFeatures(filename)
  } else {
    tdb <- makeTranscriptDbFromBiomart(biomart=biomart, dataset=dataset)
    dir.create(path=cache.dir, recursive=TRUE, showWarnings=FALSE)
    saveFeatures(tdb, file=filename)
  }
  tdb
}

## ## For SPP package
## read.bed.tags <- function(bedfile) {
##   f <- file(bedfile)
##   df <- sqldf("select * from f", dbname=tempfile(), file.format=list(header=FALSE, row.names=FALSE, sep="\t"))
##   close(f)

##   seqnames <- factor(df[[1]])
##   strand <- df[[6]]
##   pos1 <- as.integer(df[[2]])
##   pos2 <- as.integer(df[[3]])
##   score <- as.integer(df[[5]])

##   tag.locations <- ifelse(strand == "+", pos1, pos2)
##   list(tags=split(tag.locations, seqnames),
##        quality=split(score, seqnames))
## }

## read.tags <- function(file, type=NULL, ...) {
##   if (is.null(type)) {
##     type <- file.extension(file)
##   }
##   type <- type[[1]]
##   tag.reader <- get(sprintf("read.%s.tags", type))
##   tag.reader(file, ...)
## }

do.cached <- function(expr, cache.file) {
  if (file.exists(cache.file)) {
    eval(readRDS(cache.file))
  } else {
    result <- expr
    saveRDS(result, file=cache.file)
    result
  }
}

`%=%` <- function(lhs, rhs) {
  lhs <- as.character(as.symbol(substitute(lhs)))
  cache.file <- sprintf(".%s.RDS", lhs)
  delayedAssign(lhs, do.cached(rhs, cache.file), assign.env=parent.frame(1))
}

## Alternate version that preserves names
mlply <- function (.data, .fun = NULL, ..., .expand = TRUE, .progress = "none",
                   .parallel = FALSE)
{
  if (is.matrix(.data) & !is.list(.data))
    .data <- plyr:::.matrix_to_df(.data)
  f <- splat(.fun)
  result <- alply(.data = .data, .margins = 1, .fun = f, ..., .expand = .expand,
                  .progress = .progress, .parallel = .parallel)
  names(result) <- rownames(.data)
  result
}

names.to.options <- function(names, dashes=TRUE, dot.to.dash=!as.is, as.is=FALSE) {
  ## Sanity check
  stopifnot(all(str_detect(names, "[^-.]")))
  stopifnot(all(str_length(names) > 0))

  if (dot.to.dash) {
    option.strings <- str_replace(names, "\\.", "-")
  } else {
    option.strings <- names
  }
  if (identical(dashes, TRUE)) {
    ## dashes=TRUE (default) means choose automatically.
    ## Two dashes for multi-char options and one for single-char
    dashes <- ifelse(str_length(option.strings) > 1, 2, 1)
  } else if (!dashes) {
    ## dashes=FALSE means do not prepend dashes at all (equivalent to dashes=0)
    dashes <- 0
  } else {
    ## Anything else is an integer vector
    dashes <- as.integer(as.numeric(dashes))
    if (!all(dashes %in% 0:2)) {
      warning("Weird number of dashes requested")
    }
  }
  dash.strings <- str_dup("-", dashes)
  str_c(dash.strings, option.strings)
}

get.option.vector <- function(args, ...) {
  stopifnot(all(Map(length, args) == 1))
  if (is.null(names(args))) {
    names(args) <- rep(NA, length(args))
  }
  foreach(name=(if (is.null(names(args))) "" else names(args)), arg=args, .combine=c) %do% {
    if (name == "" || is.na(name)) {
      ## Positional argument
      stopifnot(!is.logical(arg))
      as.character(arg)
    } else {
      ## Option
      option.string <- names.to.options(name, ...)
      if (is.logical(arg)) {
        ## Boolean flag
        if (arg) {
          option.string
        }
      } else {
        ## Option with string argument
        c(option.string, as.character(arg))
      }
    }
  }
}

command.to.string <- function(command, args) {
  atoms <- c(command, args)
  str_c(unlist(llply(atoms, shQuote)), collapse=" ")
}

command.exists <- function(command) {
  system2("which", command, stdout=FALSE) == 0
}

do.in.directory <- function(dir, expr) {
  oldpwd <- if (!is.null(dir)) setwd(dir) else NULL
  tryCatch(expr, finally=setwd(oldpwd))
}

.abspath.single <- function(..., fsep = .Platform$file.sep) {
  ## Eliminate everything before the last path component that is
  ## absolute, if any
  pathcomps <- c(...)
  last.abs.pos <- max(which(str_sub(pathcomps, 1, 1) == fsep), -Inf)
  if (last.abs.pos == -Inf) {
    ## No absolute path components, so expand from pwd
    do.call(file.path, as.list(c(getwd(), pathcomps)))
  } else {
    do.call(filePath, as.list(pathcomps[last.abs.pos:length(pathcomps)]))
  }
}

## Like file.path, but cuts off everything before a root path, like
## Python's os.path.abspath.
abspath <- function(..., fsep = .Platform$file.sep) {
  if (.Platform$OS.type != "unix") {
    warning("abspath function is only written for UNIX paths. This will probably fail.")
  }
  unlist(Map(.abspath.single, ..., fsep=fsep))
}

make.command.runner <- function(default.path, use.intemp.by.default=FALSE) {
  function(args=list(), ..., dir=NULL, use.intemp=use.intemp.by.default, path=default.path) {
    stopifnot(command.exists(path))
    cmd <- path
    args <- get.option.vector(args)
    if (use.intemp) {
      if (!command.exists("intemp")) {
        warning("Intemp requested but not available. Continuing without it.")
      } else {
        ##     intemp --temp_dir . --preserve_temp_dir never --overwrite --
        args <- c(get.option.vector(list(temp_dir=".",
                                         preserve_temp_dir="never",
                                         overwrite=TRUE,
                                         "--")),
                  cmd, args)
        cmd <- "intemp"
      }
    }
    if (is.null(dir)) {
      dir <- getwd()
    }
    do.in.directory(dir, {
      message(sprintf("Command: %s", command.to.string(cmd, args)))
      system2(cmd, args=shQuote(args), ...)
    })

  }
}

split.ranged.data <- function(rd, chunk.size) {
  num.chunks <- ceiling(nrow(rd) / chunk.size)
  chunk.vector <- rep(seq(num.chunks), each=chunk.size, length.out=nrow(rd))
  split(rd, chunk.vector)
}

annotatePeakInBatch.parallel <- function(myPeakList, ..., chunkSize=1000) {
  tsmsgf("Splitting peak list into chunks of size %s", chunkSize)
  chunked.peakList <- split.ranged.data(myPeakList, chunkSize)
  tsmsgf("Split peak list into %s chunks", length(chunked.peakList))
  foreach(chunk.id=seq_along(chunked.peakList),mpl=chunked.peakList, .combine=rbind) %dopar% {
    tsmsgf("Processing chunk %s", chunk.id)
    annotatePeakInBatch(mpl, ...)
  }
}

strip.lcPrefix <- function(s) {
  str_sub(s, start=str_length(lcPrefix(s)) + 1)
}

## Return the highest-numbered genome provider version. For example,
## "mm10" > "mm9".
max.provider.version <- function(versions) {
  if (is.numeric(versions)) {
    return(max(versions))
  }
  v.numbers <- suppressWarnings(as.numeric(strip.lcPrefix(versions)))
  if (any(is.na(v.numbers))) {
    ## Differing suffixes are non-numeric, so choose max
    ## lexicographically
    versions[order(versions, decreasing=TRUE)[1]]
  } else {
    ## Differing suffixes are numeric, as in "mm9" vs "mm10", so
    ## choose max numerically.
    versions[order(v.numbers, decreasing=TRUE)[1]]
  }
}

get.genome <- function(genome) {
  ## First try as a package name
  ag <- available.genomes(TRUE)
  lcgenome <- tolower(genome)
  if (genome %in% ag$pkgname) {
    if (library(genome, character.only=TRUE, logical.return=TRUE)) {
      x <- sprintf("package:%s", genome)
      get(ls(pos=x)[[1]], pos=x)
    } else {
      stop("You need to install the ", genome, " package to use this genome\nTry: biocLite(\"", genome, "\")")
    }
  } else if (sprintf("BSgenome.%s", genome) %in% ag$pkgname) {
    get.genome(sprintf("BSgenome.%s", genome))
  } else if (lcgenome %in% tolower(ag$provider_version)) {
    ## Figure out the corresponding package name
    pkg <- ag$pkgname[lcgenome == tolower(ag$provider_version)][[1]]
    get.genome(pkg)
  } else if (lcgenome %in% tolower(ag$organism)) {
    ## Get max provider version for that organism
    pversion <- max.provider.version(ag$provider_version[lcgenome == tolower(ag$organism)])
    get.genome(pversion)
  } else {
    stop(sprintf("Unknown genome: %s", genome))
  }
}

group.columns <- function(df, ..., factors=list()) {
  factors.and.args <- c(as.list(factors), list(...), lex.order=TRUE, drop=TRUE)
  Map(data.frame, split(as.list(df), do.call(interaction, factors.and.args)), row.names=list(rownames(df)), check.names=FALSE)
}

group.and.pair <- function(df, grouping.factors=list(TRUE), pairing.factor=1:2, select=TRUE) {
  if (is.list(select)) {
    select <- Reduce(`&`, select)
  }
  select <- as.logical(select)

  if (!is.list(grouping.factors)) {
    grouping.factors <- list(grouping.factors)
  }

  grouping <- do.call(interaction, c(grouping.factors, lex.order=TRUE))

  pairing.factor <- pairing.factor[select]
  grouping <- grouping[select]
  df <- df[select]

  stopifnot(length(unique(pairing.factor)) == 2)

  pairs <- split(as.list(df), grouping)
  pair.labels <- rep(split(pairing.factor, grouping), length.out=length(pairs))
  for (i in 1:length(pairs)) {
    stopifnot(length(pairs[[i]]) == 2)
    stopifnot(length(unique(pair.labels[[i]])) == 2)
    names(pairs[[i]]) <- as.character(pair.labels[[i]])
  }

  pairs
}

#function(counts=counts, mymain="Venn Diagram", mysub="default", setlabels="default", yoffset=seq(0,10,by=0.34), ccol=rep(1,31), colmode=1, lcol=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), lines=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), mylwd=3, diacol=1, type="ellipse", ccex=1.0, lcex=1.0, sepsplit="_", ...) {

source("overLapper.R")

## Takes a data frame of logical columns
myVennPlot <- function(df, ...) {
  df <- as.data.frame(df)
  setlist <- llply(df, function(x) row.names(df)[x])
  ollist <- overLapper(setlist, type="vennsets")
  counts <- sapply(ollist$Venn_List, length)
  arglist <- rename(list(...), c("main"="mymain", "sub"="mysub"))
  arglist$counts <- counts
  if (!"mysub" %in% names(arglist)) {
    arglist$mysub <-   str_c(sprintf("Total counts: All: %s; ", length(unique(unlist(setlist)))),
                             str_c(sprintf("S%s: %s", seq_along(setlist), sapply(setlist, length)),
                                   collapse="; "))
  }
  do.call(vennPlot, arglist)
}

makeBrowserSession <- function(tracks, genome="hg19", browser="UCSC") {
  tsmsgf("Creating session for %s browser", browser)
  session <- browserSession(browser)
  foreach(tname=names(tracks)) %do% {
    tsmsgf("Adding track %s to session", tname)
    track(session, tname) <- tracks[[tname]]
  }
  session
}

withCores <- function(n, expr) {
  prev.cores <- options(cores=min(n, parallel:::detectCores()))$cores
  on.exit(options(cores=prev.cores))
  expr
}

delete.zero.columns <- function(x) {
  x[,colSums(abs(x)) > 0]
}

reorder.columns <- function(df, start=character(0), end=character(0)) {
  if (is.character(start)) {
    start <- na.omit(match(start, names(df)))
  } else if (is.logical(start)) {
    start <- which(start)
  } else {
    start <- as.numeric(start)
  }
  if (is.character(end)) {
    end <- na.omit(match(end, names(df)))
  } else if (is.logical(end)) {
    end <- which(end)
  } else {
    end <- as.numeric(end)
  }
  if (any(start %in% end)) {
    stop("Cannot move the same column to the start and end.")
  }
  mid <- setdiff(1:ncol(df), c(start, end))
  df[c(start, mid, end)]
}

## perldoc -f quotemeta
quotemeta <- function (string)
{
  str_replace_all(string, "(\\W)", "\\\\\\1")
}

str_stripprefix <- function(string, prefix) {
  str_replace(string, str_c("^", quotemeta(prefix)), "")
}

retrieveParallelResult <- function(p) {
  result <- tryCatch(mccollect(p)[[1]], error=function(...) NULL)
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

`%=%` <- function(lhs, rhs) {
  lhs <- as.character(as.symbol(substitute(lhs)))
  rhs.parallel <- mcparallel(rhs)
  delayedAssign(lhs, retrieveParallelResult(rhs.parallel), assign.env=parent.frame(1))
}

## Alternate version that preserves names
mlply <- function (.data, .fun = NULL, ..., .expand = TRUE, .progress = "none",
                   .parallel = FALSE)
{
  if (is.matrix(.data) & !is.list(.data))
    .data <- plyr:::.matrix_to_df(.data)
  f <- splat(.fun)
  result <- alply(.data = .data, .margins = 1, .fun = f, ..., .expand = .expand,
                  .progress = .progress, .parallel = .parallel)
  names(result) <- rownames(.data)
  result
}

table.to.granges <- function(table,
                             seqnames.column="seqnames", start.column="start",
                             end.column="end", strand.column="strand",
                             start.offset=1, end.offset=0, seqlengths=NULL) {
  special.columns <- c(seqnames.column, start.column, end.column, strand.column)
  stopifnot(all(special.columns %in% names(table)))
  metadata.columns <- setdiff(names(table), special.columns)

  if (is.character(seqlengths) && length(seqlengths) == 1) {
    seqlengths <- seqlengthsForUCSCGenome(seqlengths)
  }

  do.call("GRanges",
          c(list(seqnames=Rle(as.factor(table[[seqnames.column]])),
                 ranges=
                 IRanges(start=as.integer(table[[start.column]]+start.offset),
                         end=as.integer(table[[end.column]]+end.offset)),
                 strand=Rle(as.factor(table[[strand.column]]))),
            if (!is.null(seqlengths)) list(seqlengths=seqlengths) else list(),
            Map(Rle, table[metadata.columns])))
}

granges.to.dataframe <- function(gr, ignore.strand=FALSE, include.width=FALSE) {
  df.columns <- list(row.names=names(gr),
                     chr=seqnames(gr),
                     start=start(gr),
                     end=end(gr))
  if (include.width) {
    df.columns$width <- df.columns$end - df.columns$start + 1
  }
  if (!ignore.strand) {
    df.columns$strand <- strand(gr)
  }
  df.columns <- c(df.columns, as.list(elementMetadata(gr)))
  do.call(DataFrame, df.columns)
}

merge.nearby.ranges <- function(x, merge.distance) {
  x <- extend.ranges(x, merge.distance, side="end")
  x <- reduce(x)
  x <- extend.ranges(x, -merge.distance, side="end")
  x
}

extend.ranges <- function(x, amount, side="both") {
  stopifnot(side %in% c("start", "end", "both"))
  if (side %in% c("start", "both")) {
    start(x) <- start(x) - amount
  }
  if (side %in% c("end", "both")) {
    end(x) <- end(x) + amount
  }
  x
}

mc.cleanup <- function() {
  parallel:::mccollect(wait=FALSE)
  all <- parallel:::children()
  if (length(all)) {
    parallel:::mckill(all)
    parallel:::mccollect(all)
  }
}

## Returns "+" for same strand, "-" for opposite strand, "*" if either
## input is "*"
compare.strand <- function(s1, s2) {
  unstranded <- s1 == "*" | s2 == "*"
  same.strand <- s1 == s2
  relstrand <- ifelse(unstranded, "*",
                      ifelse(same.strand, "+", "-"))
}

duplicate.count <- function(v) {
  as.vector(table(v)[v])
}

## Returns a copy of GRanges object with all strand set to "*"
unstranded <- function(gr) {
  strand(gr) <- "*"
  gr
}

textConnectionFromLines <- function(lines, linesep="\n") {
  textConnection(str_c(sprintf("%s%s", lines, linesep), collapse=""))
}
generate.null.ranges <- function(y) {
  x <- GRanges(seqnames=names(seqlengths(y)), ranges=IRanges(1,1), strand="*", seqlengths=seqlengths(y))
  for (i in names(elementMetadata(y))) {
    elementMetadata(x)[[i]] <- as(NA, class(as.vector(elementMetadata(y)[[i]])))
  }
  x
}

select.nearest <- function(x, y) {
  y <- append(y, generate.null.ranges(y))
  y[nearest(x,y)]
}

annotate.by.granges <- function(peaks, gr, annot.columns) {
  for (i in names(elementMetadata(gr))) {
    elementMetadata(gr)[[i]] <- as.vector(elementMetadata(gr)[[i]])
  }
  nearest.ranges <- select.nearest(peaks, gr)
  nearest.distance <- distance(peaks, nearest.ranges, ignore.strand=TRUE)
  in.ranges <- Rle(nearest.distance == 0)
  annot.data <- elementMetadata(nearest.ranges)[annot.columns]
  if (!is.null(names(annot.columns))) {
    names(annot.data) <- names(annot.columns)
  }
  DataFrame(overlap=in.ranges, distance=nearest.distance, annot.data)
}

getCenters <- function(gr) {
  start(gr) <- end(gr) <- floor((start(gr) + end(gr)) / 2)
  gr
}

make.locus.names <- function(gr) {
  sprintf("%s:%s-%s(%s)", seqnames(gr), start(gr), end(gr), strand(gr))
}

slurm.expand.nodelist <- function(nl) {
  if (nl == "") {
    return("")
  }
  commasplit <- function(s) str_split(s, ",")
  expand.bracket.expression <- function(expr) {
    if (expr == "") {
      return("")
    }
    subexprs <- unlist(commasplit(expr))
    range.matches <- str_match(subexprs, "^([^-]+)-([^-]+)")
    foreach(expr=subexprs, start=range.matches[,2], end=range.matches[,3], .combine=c) %do% {
      if (is.na(start)) {
        expr
      } else {
        as.character(seq(as.numeric(start), as.numeric(end)))
      }
    }
  }
  bracketsplit <- function(expr) {
    str_match(expr, "^([^[]*)(?:\\[([^]]+)\\](.*))?$")[1,-1]
  }
  bsplit <- bracketsplit(nl)
  prefixes <- str_c(bsplit[1], expand.bracket.expression(bsplit[2]))
  suffixes <- slurm.expand.nodelist(bsplit[3])
  str_c(rep(prefixes, each=length(suffixes)),
        rep(suffixes, times=length(prefixes)))

}

get.slurm.hostlist <- function() {
  nodelist <- slurm.expand.nodelist(Sys.getenv("SLURM_NODELIST"))
  if (all(nodelist == "")) {
    stop("Invalid SLURM_NODELIST environment variable. Is this script running under SLURM?")
  }
  cpus.per.node <- as.numeric(str_match(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"), "^[0-9]+")[1,1])
  if (is.na(cpus.per.node)) {
    stop("Invalid SLURM_JOB_CPUS_PER_NODE environment variable. Is this script running under SLURM?")
  }
  rep(nodelist, cpus.per.node)
}

makeNiceCluster.byslurm <- function(nice=TRUE) {
  cluster.host.list <- get.slurm.hostlist()
  (if (nice) makeNiceCluster else makeCluster)(cluster.host.list)
}

merge.defaults.into.list <- function(x, defaults) {
  combined <- c(x, defaults)
  combined[!duplicated(names(combined))]
}

curry <- function(f, args) {
  args <- as.list(args)
  if (length(args) == 0) {
    f
  } else {
    function(...) do.call(f, c(list(...), as.list(args)))
  }
}

## Return what print(...) would print, as a single character string
sprint <- function(...) {
  str_c(capture.output(print(...)),
        collapse="\n")
}

## If llply is called on an unnamed character vector, use the
## character vector itself for the names.
llply <- function(.data, ...) {
  if (is.character(.data) && is.null(names(.data))) {
    names(.data) <- .data
  }
  plyr::llply(.data, ...)
}

get.dots.uneval <- function() as.list(substitute((...), env = parent.frame()))[-1]

eval.par <- function(..., llply.args=list()) {
  expressions <- get.dots.uneval()
  do.call(llply, c(list(expressions, eval, .parallel=TRUE),
                   llply.args))
}

assign.par <- function(..., llply.args=list()) {
  expressions <- get.dots.uneval()
  ## Filter out unnamed expressions
  expressions <- expressions[names(expressions) != ""]
  results <- do.call(llply, c(list(expressions, eval, .parallel=TRUE),
                              llply.args))
  for (varname in names(results)) {
    assign(varname, results[[varname]], env=parent.frame())
  }
}

assignFrom <- function(env, verbose=TRUE) {
    stopifnot(!is.null(names(env)))
    for (n in names(env)) {
        if (n != "") {
            if (verbose)
                message("Assigning ", n)
            assign(n, env[[n]], env=parent.frame())
        }
    }
}

## Convert CharacterList to character
clist.to.character <- function(ccl, sep=",") {
  setNames(unlist(lapply(ccl, function(ch) str_c(na.omit(ch), collapse=sep))),
           names(ccl))
}

get.gene.lengths <- function(txdb, summary.fun=max) {
  ## Get the length of each transcript
  tx.ex <- exonsBy(txdb, "tx")
  tx.len <- sapply(width(tx.ex), sum)
  ## Group transcript lengths by gene
  gn.tx <- transcriptsBy(txdb, "gene")
  gn.lens <- split(tx.len[as.character(mcols(gn.tx@unlistData)$tx_id)],
                   rep(names(gn.tx), elementLengths(gn.tx)))
  ## Summarize lengths for each gene
  sapply(gn.lens, summary.fun)
}

gc.content <- function(dna) {
  afreq <- alphabetFrequency(dna)
  setNames(rowSums(afreq[,c("G", "C")] / rowSums(afreq)),
           names(dna))
}

get.gene.gc <- function(txdb, genome, summary.fun=mean) {
  ## Get the sequence content of each transcript
  tx.seq <- extractTranscriptsFromGenome(genome, txdb)
  tx.gc <- gc.content(tx.seq)

  ## Group transcript GCs by gene
  gn.tx <- transcriptsBy(txdb, "gene")
  gn.gcs <- split(tx.gc[as.character(mcols(gn.tx@unlistData)$tx_name)],
                   rep(names(gn.tx), elementLengths(gn.tx)))
  ## Summarize lengths for each gene
  sapply(gn.gcs, summary.fun)
}

## get.gene.lengths.from.gff <- function(gff, summary.fun=max) {
## }

modifyenv <- function(obj, newvalues) {
  newenv <- as.environment(newvalues)
  parent.env(newenv) <- environment(obj)
  environment(obj) <- newenv
  obj
}

withOptions <- function(opts, expr) {
  oldvalues <- options(as.list(opts))
  on.exit(options(as.list(oldvalues)))
  expr
}

withRecovery <- function(expr) {
  withOptions(list(error=recover),
              expr)
}

fix.model.matrix.names <- function(mm) {
  colnames(mm) <- str_replace_all(colnames(mm), fixed("(Intercept)"), "Intercept")
  colnames(mm) <- make.names(colnames(mm), unique=TRUE)
  mm
}

model.matrix.default <- function(...) {
  mm <- stats::model.matrix.default(...)
  fix.model.matrix.names(mm)
}

getlist <- function(varnames) {
  mget(varnames, envir=parent.frame())
}

geom.mean <- function(x) exp(mean(log(x)))
geom.mean.norm <- function(x) x / geom.mean(x)

export.txdb.to.gtf <- function(txdb, con) {
    alltx <- transcripts(txdb)
    txnames <- c()
    txnames[mcols(alltx)$tx_id] <- mcols(alltx)$tx_name
    names(txnames) <- as.character(mcols(alltx)$tx_id)
    tx.tss <- str_c(as.vector(seqnames(alltx)),
                    as.vector(strand(alltx)),
                    ifelse(as.vector(strand(alltx)) == "+",
                           start(alltx), end(alltx)),
                    sep=":")
    tx.tss <- factor(tx.tss)
    levels(tx.tss) <- sprintf("TSS%s", seq_along(levels(tx.tss)))
    names(tx.tss) <- as.character(mcols(alltx)$tx_id)

    txbg <- transcriptsBy(txdb, "gene")
    mcols(txbg@unlistData)$gene_id <- rep(names(txbg), elementLengths(txbg))
    txgenes <- with(mcols(unlist(txbg)), setNames(gene_id, tx_id))
    geneless.transcripts <- setdiff(as.character(mcols(alltx)$tx_id), names(txgenes))
    txgenes[geneless.transcripts] <- sprintf("NOGENE%s", seq_along(geneless.transcripts))

    ex <- exonsBy(txdb, "tx")
    mcols(ex@unlistData) <-
        within(mcols(ex@unlistData), {
            transcript_number <- rep(names(ex), elementLengths(ex))
            transcript_id <- txnames[transcript_number]
            tss_id <- tx.tss[transcript_number]
            exon_number <- exon_rank
            gene_id <- txgenes[transcript_number]
            type <- "exon"
        })
    export(unlist(unname(ex)), con, format="GTF")
}

readChrominfoFromFasta <- function(fafile, is_circular=FALSE) {
    fa <- readFasta(fafile)
    data.frame(chrom=as.character(id(fa)),
               length=width(fa),
               is_circular=is_circular)
}

dropEmptyCols <- function(df) {
    isNonEmpty <- function(col) {
        if (is(col, "CompressedList")) {
            length(unlist(col)) > 0
        } else {
            ! all(is.na(col))
        }
    }
    df[unlist(Map(isNonEmpty, df))]
}

callf <- function(func, place, ...) {
    args <- as.list(match.call())[-1]
    inner.call <- as.call(args)
    names(inner.call)[2] <- NA
    assign.call <- as.call(list(quote(`<-`),
                                args$place,
                                as.call(inner.call)))
    print(assign.call)
    eval(assign.call, envir=parent.frame())
}
