#!/bin/bash 
#SBATCH --job-name=sstacks
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
# run `sstacks`
############################

# sstacks is the third step of the stacks de novo pipeline

# load software--------------------------------------------------------------------------------
module load stacks/2.53

# input, output files, directories-------------------------------------------------------------
INDIR=../results/stacks/denovo
POPMAP=../meta/popmap_total.txt

# run sstacks----------------------------------------------------------------------------------
sstacks -P $INDIR -M $POPMAP -p 10

date