#!/bin/bash
# -*- mode:R -*-
#PBS
#PBS -j oe -o idr-02-downsample-bam.log
module load R/3.1.0;
cd $PBS_O_WORKDIR

Rscript --version
Rscript - <<'EOF'

library(xlsx)
library(plyr)
library(stringr)

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

fixdfchar <- function(df) {
    for (i in seq_along(df)) {
        if (is.factor(df[[i]]) && !any(duplicated(df[[i]]))) {
            df[[i]] <- as.character(df[[i]])
        }
    }
    df
}
splitargs <- function(...) {
    all.args <- unlist(list(...))
    unname(unlist(str_split(all.args, "\\s+")))
}

makecmd <- function(...) {
    args <- as.character(unlist(list(...)))
    cmd <- str_c(shQuote(args), collapse=" ")
}

make.picard.cmd <- function(jarfile, args=c(), ..., java.path="java", java.opts=c()) {
    args <- c(args, ...)
    if (is.null(names(args)) || any(is.na(names(args))) || any(names(args) == "") || length(args) == 0)
        stop("Invalid args")
    makecmd(java.path, java.opts, "-jar", jarfile, sprintf("%s=%s", names(args), args))
}

whoami <- function() {
    system2("whoami", stdout=TRUE)
}

job.name.running <- function(name) {
    unlist(
        mclapply(name,
                 function(n)
                 length(system2("qselect", c("-N", n, "-s", "EHQRTW"), stdout=TRUE)) > 0,
                 mc.cores=50))
}

submit.job <- function(cmds, qsub.opts="", wd=getwd(),
                       name=str_replace_all(cmd, "[^A-Za-z0-9_-=+:,.!@#$%^&*()[]{}<>?|]", "_"),
                       logfile,
                       modules=c(), dry.run=FALSE,
                       allow.duplicate.name=FALSE) {
    if (!allow.duplicate.name && !is.null(name) && job.name.running(name)) {
        warning("Not starting duplicate job ", name)
        return()
    }
    if (!missing(logfile)) {
        qsub.opts <- c(qsub.opts, str_c("-j oe -o ", logfile))
    }
    temp <- tempfile()
    text <- c("#!/bin/sh",
              sprintf("#PBS -w %s -N %s", wd, name),
              str_c("#PBS ", str_c(qsub.opts, collapse=" ")),
              str_c("cd $PBS_O_WORKDIR"),
              sprintf("module load %s", modules),
              cmds)
    if (dry.run) {
        writeLines(text)
    }
    else {
        writeLines(text, con=temp)
        system2("qsub", args=temp)
        file.remove(temp)
    }
}

string.ends.with <- function(str, suffix) {
    str_sub(str, -str_length(suffix)) == suffix
}

file.name.as.directory <- function(f) {
    d <- file.path(f)
    if (! string.ends.with(d, .Platform$file.sep)) {
        d <- str_c(d, .Platform$file.sep)
    }
    d
}

wd <- "/gpfs/group/salomon/userdata/ryan/Projects/sarah-cd4"
setwd(wd)

mergedir <- file.path(wd, "merged_bam_files")
sample.prob <- 0.1

tsmsg("Assembling downsample targets")

stypes <- c("input", "H3K4me2", "H3K4me3", "H3K27me3")
jobnames <- setNames(sprintf("DS:%s", stypes), stypes)

big.bam.files <- setNames(sprintf("%s/%s.bam", mergedir, stypes), stypes)
downsampled.bam.files <- setNames(sprintf("%s/%s_downsample.bam", mergedir, stypes), stypes)
logfiles <- setNames(sprintf("%s/%s.log", mergedir, jobnames), stypes)

tsmsg("Generating commands")

downsamplesam.jarfile <- Sys.which("DownsampleSam.jar")[1]
stopifnot(length(downsamplesam.jarfile) && !is.na(downsamplesam.jarfile) && downsamplesam.jarfile != "")

make.downsample.command <- function(input, output) {
    make.picard.cmd(jarfile=downsamplesam.jarfile, I=input, O=output, P=sample.prob)
}

downsample.cmds <- mapply(make.downsample.command, input=big.bam.files, output=downsampled.bam.files)

for (i in stypes) {
    of <- downsampled.bam.files[[i]]
    ofq <- shQuote(of)
    wdq <- shQuote(wd)
    jn <- jobnames[i]
    dsc <- downsample.cmds[i]
    cmds <- c(
        sprintf("cd %s", shQuote(wd)),
        sprintf("echo \"Starting DownsampleSam run %s\"", jn),
        dsc,
        sprintf("
if [ $? != 0 ]; then
  echo \"DownsampleSam run %s exited with non-zero return code. Deleting output file.\"
  rm -f %s
fi",
                jn, ofq),
        sprintf("
if [ -e %s ]; then
  echo \"Finished downsampling. Now indexing\"
  bamtools index -in %s;
  echo \"DownsampleSam run %s completed successfully on `date`\"
else
  echo \"DownsampleSam run %s FAILED on `date`\"
  rm -f %s
fi",
                ofq, ofq, jn, jn, ofq))
    if (job.name.running(jn)) {
        tsmsg("Skipping job ", jn, " because it is already currently running.")
        next()
    } else if (!file.exists(big.bam.files[i])) {
        tsmsg("Skipping job ", jn, " because of missing input files.")
        next()
    } else if (file.exists(of)) {
        tsmsg("Skipping job ", jn, " because it is already complete.")
        next()
    }
    tsmsg("Submitting job ", jn)
    submit.job(cmds, modules=c("samtools", "bamtools", "picard/1.92"),
               name=jn,
               logfile=logfiles[i], wd=wd,
               dry.run=FALSE)
}

EOF <- NULL
EOF
