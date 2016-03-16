#!/bin/bash
# -*- mode:R -*-
#PBS
#PBS -j oe -o idr-03-peakcall.log
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

qsub.opts <- "-l mem=15gb,walltime=48:00:00"

tsmsg("Loading sample data")

expdata <- read.xlsx("sampledata.xlsx", 1)
expdata <- fixdfchar(expdata)
rownames(expdata) <- expdata$Sample
expdata$Lane <- factor(expdata$Lane)
expdata$Timepoint <- factor(expdata$Timepoint, levels=str_c("D", c(0, 1, 5, 14)))
expdata$Celltype <- factor(expdata$Celltype, levels=c("Naive", "Memory"))
expdata$Sampletype <- factor(expdata$Sampletype, levels=c("input", "H3K4me2", "H3K4me3", "H3K27me3"))

bamdir <- file.path(wd, "bam_files")
expdata$bampath <- file.path(bamdir, expdata$BAM)
stopifnot(all(file.exists(expdata$bampath)))


mergedir <- file.path(wd, "merged_bam_files")

tsmsg("Assembling non-input merge targets")

expdata.noninput <- droplevels(expdata[expdata$Sampletype != "input",])
rownames(expdata.noninput) <- with(expdata.noninput, make.names(Sampletype:Celltype:Timepoint:Donor))
groups.by.stype <- with(expdata.noninput, split(bampath, Sampletype))
names(groups.by.stype) <- str_c(names(groups.by.stype))
groups.by.cond <- with(expdata.noninput, split(bampath, Sampletype:Celltype:Timepoint))
groups.by.stype.and.donor <- with(expdata.noninput, split(bampath, Sampletype:Donor))

merge.targets <- c(groups.by.stype, groups.by.stype.and.donor, groups.by.cond)
merge.targets <- setNames(llply(merge.targets, as.character), make.names(names(merge.targets)))
mergenames <- setNames(names(merge.targets), names(merge.targets))

single.chip.bams <- setNames(expdata.noninput$bampath, rownames(expdata.noninput))
merged.input.bam <- file.path(wd, "merged_bam_files", "input_downsample.bam")
stopifnot(file.exists(merged.input.bam))
merged.chip.bams <- setNames(sprintf("%s/%s.bam", mergedir, mergenames), mergenames)
## If a downsampled version exists, use it
downsampled.bams <- setNames(sprintf("%s/%s_downsample.bam", mergedir, mergenames), mergenames)
downsampled <- file.exists(downsampled.bams)
merged.chip.bams[downsampled] <- downsampled.bams[downsampled]

stopifnot(all(file.exists(c(merged.input.bam, merged.chip.bams))))

tsmsg("Generating commands")

peakcall.dir <- file.path(wd, "peak_calls")
dir.create(peakcall.dir, recursive=TRUE, showWarnings=FALSE)

peakcall.jobs <- data.frame(chip.bam=c(single.chip.bams, merged.chip.bams), stringsAsFactors=FALSE)
peakcall.jobs <- within(peakcall.jobs, {
    jobname <- sprintf("PC:%s", row.names(peakcall.jobs))
    outbase <- file.path(peakcall.dir, row.names(peakcall.jobs))
    logfile <- file.path(peakcall.dir, sprintf("%s.log", jobname))
    outfile <- sprintf("%s_peaks.encodePeak", outbase)
    cmd <- with(peakcall.jobs, {
        mapply(makecmd,
               "macs2", "callpeak",
               "--treatment", chip.bam,
               "--control", merged.input.bam,
               "--name", outbase,
               "--gsize", "hs",
               "--pvalue", "0.1",
               ## Broad peak cutoff needs to be bigger than main
               ## p-value cutoff
               "--broad-cutoff", "0.2",
               ## Auto-determination of shift size isn't working, so
               ## we just set it to ~nucleosome size
               "--nomodel", "--shiftsize", "74",
               "--bdg",
               "--to-large",
               "--broad",
               "--call-summits")
    })
})

for (i in rownames(peakcall.jobs)) {
    attach(peakcall.jobs[i,])
    ofq <- shQuote(outfile)
    cmds <- c(
        sprintf("cd %s", shQuote(wd)),
        sprintf("mkdir -p %s", shQuote(peakcall.dir)),
        sprintf("echo \"Starting MACS peak call run %s\"", jobname),
        cmd,
        sprintf("
if [ $? != 0 ]; then
  echo \"MACS peak call run %s FAILED: MACS exited with non-zero return code. Deleting output files.\"
  rm -f %s.*
", jobname, shQuote(outbase)),
        sprintf("
elif [ -e %s ]; then
  echo \"MACS peak call run %s completed successfully on `date`\"
else
  echo \"MACS peak call run %s FAILED to produce output on `date`\"
fi", ofq, jobname, jobname))
    if (job.name.running(jobname)) {
        tsmsg("Skipping job ", jobname, " because it is already currently running.")
        detach(); next()
    } else if (!all(file.exists(c(chip.bam, merged.input.bam)))) {
        tsmsg("Skipping job ", jobname, " because of missing input files.")
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
               dry.run=FALSE)
    detach();
}

EOF <- NULL
EOF
