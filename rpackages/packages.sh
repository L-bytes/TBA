#!/bin/bash
#SBATCH -N 1
#SBATCH -t 48:00:00

module load 2022
module load R

count=ls $HOME/rpackages/*.tar.gz | wc -l
directory= ls -d */ | wc -l

while [ $directory -lt $count ]
do
  for filename in $HOME/rpackages/*.tar.gz
  do
    export R_LIBS=$HOME/rpackages:$R_LIBS
    R CMD INSTALL -l $HOME/rpackages "$filename"
  done
done