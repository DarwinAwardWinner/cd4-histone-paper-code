#!/bin/bash
# -*- mode:R -*-
#PBS
#PBS -j oe -o idr-05-bca-analysis.log
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

qsub.opts <- "-l walltime=48:00:00"

peakcall.dir <- file.path(wd, "peak_calls")

idr.dir <- file.path(wd, "idr_analysis")
dir.create(idr.dir, recursive=TRUE, showWarnings=FALSE)

bca.script <- "/gpfs/home/rcthomps/src/idrCode/batch-consistency-analysis.r"

tsmsg("Loading sample data")

expdata <- read.xlsx("sampledata.xlsx", 1)
expdata <- fixdfchar(expdata)
expdata$Lane <- factor(expdata$Lane)
expdata$Timepoint <- factor(expdata$Timepoint, levels=str_c("D", c(0, 1, 5, 14)))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
expdata$Sampletype <- factor(expdata$Sampletype, levels=c("input", "H3K4me2", "H3K4me3", "H3K27me3"))
expdata$Sample <-
    with(expdata, make.names(Sampletype:Celltype:Timepoint:Donor))
expdata$Group <- make.names(with(expdata, Sampletype:Celltype:Timepoint))

inputsamples <- expdata$Sampletype == "input"

expdata.noninput <- droplevels(expdata[!inputsamples,])

expdata.donormerged <- within(data.frame(unique(expdata.noninput[c("Sampletype", "Donor")])), {
    Sample <- make.names(Sampletype:Donor)
    Group <- Sampletype
})

expdata.allbca <- rbind(expdata.noninput[names(expdata.donormerged)], expdata.donormerged)
rownames(expdata.allbca) <- expdata.allbca$Sample
expdata.allbca$Peakfile <- file.path(peakcall.dir, sprintf("%s_peaks_sorted.encodePeak", expdata.allbca$Sample))

stopifnot(all(file.exists(expdata.allbca$Peakfile)))

idrgroups <- make.names(levels(expdata.allbca$Group))
alldonors <- levels(expdata$Donor)

tsmsg("Generating commands")

idr.jobs <- ddply(expdata.allbca, .(Group), function(df) {
    ## df <- data.frame(row.names=as.character(interaction(Group, alldonors)),
    ##                  Sample=as.character(interaction(Group, alldonors)),
    ##                  Donor=alldonors,
    ##                  Peakfile=file.path(peakcall.dir, sprintf("%s.%s_peaks_sorted.encodePeak", Group, alldonors)))
    stopifnot(all())## df <- df[file.exists(df$Peakfile),]
    stopifnot(nrow(df) >= 3)
    ## Get every pair of rows
    rowpairs <- combn(row.names(df), 2)
    r1 <- rowpairs[1,]
    r2 <- rowpairs[2,]
    data.frame(Donor1=df[r1,]$Donor,
               Donor2=df[r2,]$Donor,
               Sample1=df[r1,]$Sample,
               Sample2=df[r2,]$Sample,
               Peakfile1=df[r1,]$Peakfile,
               Peakfile2=df[r2,]$Peakfile,
               Pairname=sprintf("%s:%s.vs.%s", df$Group[1], df[r1,]$Donor, df[r2,]$Donor),
               stringsAsFactors=FALSE)
})

row.names(idr.jobs) <- idr.jobs$Pairname
idr.jobs <- within(idr.jobs, {
    jobname <- str_c("BC:", Pairname)
    logfile <- file.path(peakcall.dir, sprintf("%s.log", jobname))
    outbase <- file.path(idr.dir, make.names(Pairname))
    outfile <- str_c(outbase, "-npeaks-aboveIDR.txt")
    cmd <- mapply(makecmd, "time", "Rscript", basename(bca.script),
                  Peakfile1, Peakfile2, "-1", outbase, "0", "F", "p.value")
})

for (i in rownames(idr.jobs)) {
    attach(idr.jobs[i,])
    ofq <- shQuote(outfile)
    cmds <- c(
        sprintf("echo \"Doing BCA for job name: %s\"", jobname),
        sprintf("echo \"Started BCA at `date`\""),
        ## BCA script must be run from its own directory
        makecmd("cd", dirname(bca.script)),
        cmd,
        sprintf("
if [ $? != 0 ]; then
  echo \"BCA run %s FAILED at `date`: script exited with non-zero return code. Deleting output files.\"
  rm -f %s*
", jobname, shQuote(outbase)),
        sprintf("
elif [ -e %s ]; then
  echo \"BCA run %s completed successfully at `date`\"
else
  echo \"BCA run %s FAILED to produce output at `date`\"
  rm -f %s*
fi", ofq, jobname, jobname, shQuote(outbase)))
    if (job.name.running(jobname)) {
        tsmsg("Skipping job ", jobname, " because it is already currently running.")
        detach(); next()
    } else if (!all(file.exists(c(Peakfile1, Peakfile2)))) {
        tsmsg("Skipping job ", jobname, " because of missing input file.")
        detach(); next()
    } else if (file.exists(outfile)) {
        tsmsg("Skipping job ", jobname, " because it is already complete.")
        detach(); next()
    }
    tsmsg("Submitting job ", jobname)
    submit.job(cmds,
               qsub.opts=qsub.opts,
               name=jobname,
               logfile=logfile, wd=wd,
               dry.run=FALSE)
    detach();
}

EOF <- NULL
EOF
