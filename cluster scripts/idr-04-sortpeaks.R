#!/bin/bash
# -*- mode:R -*-
#PBS
#PBS -j oe -o idr-04-sortpeaks.log
#SBATCH -c8
module load R/3.1.0;
cd $PBS_O_WORKDIR

Rscript --version
Rscript - <<'EOF'

library(annotate)
library(xlsx)
library(plyr)
options(cores=8)
options(mc.cores=8)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(Rsubread)
## Not sure which of these is the real option. Some code seems to
## use different ones.
options(cores=parallel:::detectCores())
options(mc.cores=parallel:::detectCores())
registerDoParallel()

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
        temp <- tempfile()
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

qsub.opts <- "-l walltime=48:00:00"

tsmsg("Loading sample data")

expdata <- read.xlsx("sampledata.xlsx", 1)
expdata <- fixdfchar(expdata)
rownames(expdata) <- expdata$Sample
expdata$Lane <- factor(expdata$Lane)
expdata$Timepoint <- factor(expdata$Timepoint, levels=str_c("D", c(0, 1, 5, 14)))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
expdata$Sampletype <- factor(expdata$Sampletype, levels=c("input", "H3K4me2", "H3K4me3", "H3K27me3"))

inputsamples <- expdata$Sampletype == "input"

expdata.noninput <- droplevels(expdata[!inputsamples,])
samplenames <- with(expdata.noninput,
                    c(make.names(Sampletype:Celltype:Timepoint:Donor),
                      make.names(levels(Sampletype:Celltype:Timepoint)),
                      make.names(levels(Sampletype)),
                      make.names(levels(Sampletype:Donor))))

tsmsg("Generating commands")

peakcall.dir <- file.path(wd, "peak_calls")

num.keep <- "100000"

sortpeak.jobs <- data.frame(row.names=samplenames, sample=samplenames)
sortpeak.jobs <- within(sortpeak.jobs, {
    jobname <- sprintf("SP:%s", sample)
    infile <- file.path(peakcall.dir, sprintf("%s_peaks.encodePeak", sample))
    tempfile <- file.path(peakcall.dir, sprintf("%s_peaks_tempsort.encodePeak", sample))
    outfile <- file.path(peakcall.dir, sprintf("%s_peaks_sorted.encodePeak", sample))
    logfile <- file.path(peakcall.dir, sprintf("%s.log", jobname))
    cmd <- sprintf("sort -k 8nr,8nr %s | head -n %s > %s && mv %s %s",
                   shQuote(infile), shQuote(num.keep), shQuote(tempfile),
                   shQuote(tempfile), shQuote(outfile))
})

for (i in rownames(sortpeak.jobs)) {
    attach(sortpeak.jobs[i,])
    ofq <- shQuote(outfile)
    cmds <- c(
        sprintf("cd %s", shQuote(wd)),
        sprintf("echo \"Sorting peaks; job name: %s\"", jobname),
        sprintf("echo \"Started sorting at `date`\""),
        cmd,
        sprintf("
if [ $? != 0 ]; then
  echo \"Peak sorting run %s FAILED at `date`: sort exited with non-zero return code. Deleting output file.\"
  rm -f %s %s
", jobname, shQuote(outfile), shQuote(tempfile)),
        sprintf("
elif [ -e %s ]; then
  echo \"Peak sort run %s completed successfully at `date`\"
else
  echo \"MACS peak call run %s FAILED to produce output at `date`\"
  rm -f %s
fi", ofq, jobname, jobname, shQuote(tempfile)))
    if (job.name.running(jobname)) {
        tsmsg("Skipping job ", jobname, " because it is already currently running.")
        detach(); next()
    } else if (!file.exists(infile)) {
        tsmsg("Skipping job ", jobname, " because of missing input file.")
        detach(); next()
    } else if (file.exists(outfile)) {
        tsmsg("Skipping job ", jobname, " because it is already complete.")
        detach(); next()
    }
    tsmsg("Submitting job ", jobname)
    submit.job(cmds, modules="macs/2.0.10",
               qsub.opts=qsub.opts,
               name=jobname,
               logfile=logfile, wd=wd,
               dry.run=TRUE)
    detach();
}

EOF <- NULL
EOF
