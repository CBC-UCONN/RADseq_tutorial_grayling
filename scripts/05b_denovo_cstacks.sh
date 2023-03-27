#!/bin/bash 
#SBATCH --job-name=cstacks
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

############################
# run `cstacks`
############################

# cstacks is the second step of the stacks de novo pipeline

# load software------------------------------------------------------------
module load stacks/2.53
module load R/3.6.3

# input, output files, directories-----------------------------------------
INDIR=../results/stacks/denovo

# first create the popmap file
	# this R script creates the popmap files required for the next 
	# steps of stacks. we're using an R script so we can subset 
	# samples for cstacks. using all samples would take too long with the total 
	# dataset. the next line of code executes an R script that 
	# does this. it's commented out for the tutorial because it 
	# requires the tidyverse packages to be installed, which users 
	# may not have, but is included in the repository to completely 
	# document steps required for the analysis. 
	 
# Rscript make_popmaps.R

POPMAP=../meta/popmap_cstacks.txt

# run cstacks---------------------------------------------------------------
cstacks \
-P $INDIR \
-M $POPMAP \
-n 2 \
-p 10 