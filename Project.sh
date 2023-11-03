#!/bin/bash
#Set job requirements
#SBATCH -N 1
#SBATCH -t 27:00:00
 
#Loading modules
module load 2022
module load R
 
#Create output directory on scratch
mkdir "$TMPDIR"/"$SLURM_JOBID"

#Copy input file to scratch
cp -r "$PWD"/data "$TMPDIR"/"$SLURM_JOBID"
cp -r "$PWD"/src "$TMPDIR"/"$SLURM_JOBID"

export R_LIBS=$HOME/rpackages:$R_LIBS
 
#Execute a Python program located in $HOME, that takes an input file and output directory as arguments.
Rscript --no-save --slave "$TMPDIR"/"$SLURM_JOBID"/src/Main.R "$TMPDIR/$SLURM_JOBID" "KRAS" "$HOME" ""
 
#Copy output directory from scratch to home
mkdir "$PWD"/"$SLURM_JOBID"
cp -r "$TMPDIR"/"$SLURM_JOBID"/output "$PWD/$SLURM_JOBID"
#cp -r "$TMPDIR"/"$SLURM_JOBID"/figures "$PWD/$SLURM_JOBID/output"
#cp -r "$TMPDIR"/"$SLURM_JOBID"/rdata "$PWD/$SLURM_JOBID/output"