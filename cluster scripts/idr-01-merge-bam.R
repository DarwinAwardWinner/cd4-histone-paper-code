#!/bin/bash
# -*- mode:R -*-
#PBS -j oe -o idr-merge-bam.log
module load R/3.1.0;
cd $PBS_O_WORKDIR

Rscript --version
Rscript - <<'EOF'

library(annotate)
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
dir.create(mergedir, recursive=TRUE, showWarnings=FALSE)

tsmsg("Assembling merge targets")

groups.by.stype <- with(expdata, split(bampath, Sampletype))

expdata.noninput <- droplevels(expdata[expdata$Sampletype != "input",])
groups.by.cond <- with(expdata.noninput, split(bampath, Sampletype:Celltype:Timepoint))
groups.by.stype.and.donor <- with(expdata.noninput, split(bampath, Sampletype:Donor))

merge.targets <- c(groups.by.stype, groups.by.stype.and.donor, groups.by.cond)
merge.targets <- setNames(llply(merge.targets, as.character), make.names(names(merge.targets)))
jobnames <- setNames(sprintf("MRG:%s", names(merge.targets)), names(merge.targets))
output.files <- setNames(
    sprintf("%s/%s.bam", mergedir, names(merge.targets)),
    names(merge.targets))
logfiles <- setNames(
    sprintf("%s/%s.log", mergedir, jobnames),
    names(merge.targets))

tsmsg("Generating commands")

make.merge.command <- function(inputs, output) {
    makecmd("samtools", "merge", output, inputs)
}

merge.cmds <- sapply(names(merge.targets), function(i) {
    make.merge.command(merge.targets[[i]], output.files[[i]])
})

for (i in names(merge.cmds)) {
    of <- output.files[i]
    ofq <- shQuote(of)
    wdq <- shQuote(wd)
    jn <- jobnames[i]
    mc <- merge.cmds[i]
    cmds <- c(
        sprintf("cd %s", shQuote(wd)),
        sprintf("mkdir -p %s", shQuote(mergedir)),
        sprintf("echo \"Starting samtools merge run %s\"", jn),
        merge.cmds[[i]],
        sprintf("
if [ $? != 0 ]; then
  echo \"samtools merge run %s exited with non-zero return code. Deleting output file.\"
  rm -f %s
fi",
                jn, ofq),
        sprintf("
if [ -e %s ]; then
  echo \"Finished merging. Now indexing\"
  bamtools index -in %s;
  echo \"samtools merge run %s completed successfully on `date`\"
else
  echo \"samtools merge run %s FAILED on `date`\"
  rm -f %s
fi",
                ofq, ofq, jn, jn, ofq))
    if (job.name.running(jn)) {
        tsmsg("Skipping job ", jn, " because it is already currently running.")
        next()
    } else if (!all(file.exists(merge.targets[[i]]))) {
        tsmsg("Skipping job ", jn, " because of missing input files.")
        next()
    } else if (file.exists(of)) {
        tsmsg("Skipping job ", jn, " because it is already complete.")
        next()
    }
    tsmsg("Submitting job ", jn)
    submit.job(cmds, modules=c("samtools", "bamtools"),
               name=jn,
               logfile=logfiles[i], wd=wd,
               dry.run=FALSE)
}

EOF <- NULL
EOF
