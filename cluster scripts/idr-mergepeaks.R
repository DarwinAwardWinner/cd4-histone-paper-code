#!/bin/bash
# -*- mode:R -*-
#PBS
#PBS -j oe -o idr-mergepeaks.log
#SBATCH -c8
module load R/3.1.0;
cd $PBS_O_WORKDIR

Rscript --version
Rscript - <<'EOF'
