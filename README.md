This is the codebase for the Master thesis of M. van Laar.

The main branch is for local use, but it does not run well on a normal computer.

The mutation variations file is too big to be uploaded onto GitHub, so it can be downloaded and added to the data folder from: https://drive.google.com/file/d/1vWALr6PEFi3uZCPEj0Vo5CgmcyUzZlzK/view?usp=sharing

Instructions for local use:
1. Install R (4.2), Python (3.5), Ubuntu (16.04 LTS) and GNU parallel
   
2. Install neccessary R packages:
   BiocManager
   BiocParallel
   BioNet
   clusterprofiler
   countdata
   data.table
   DESeq2
   doParallel
   dplyr (exclude any functions that might conflict with other packages)
   enrichplot
   foreach
   ggplot2
   ggprism
   ggrepel
   ggraph
   ibb
   msigdbr
   "org.Hs.eg.db"
   pheatmap
   reshape2
   RColorBrewer
   stringr
   tidygraph
   tidyr
   "vsn"
   WGCNA

3. Run the code

You can modify the mutation and variant studied by adjusting the hit and variant variables at the top of the script.

The cluster branch is for use on a cluster computing system that runs on sbatch.

Instructions for Snellius:
1. Make alias for Python and R
  > cat ~/.modulerc.lua
  > module_alias("R", "R/4.2.1-foss-2022a")
  > module_alias("py", "Python/2.7.18-GCCcore-11.3.0-bare")

2. Install neccessary packages and dependencies by running the packages.sh script within the rpackages folder

3. Run code
  > sbatch --ntasks=1 --cpus-per-task=32 --mem-per-cpu=7000 Project.sh

You can modify the mutation and variant studied by adjusting the third and fifth arguments of the Project.sh file.
