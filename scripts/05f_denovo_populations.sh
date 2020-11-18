#!/bin/bash 
#SBATCH --job-name=populations
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G
#SBATCH --qos=general
#SBATCH --partition=general


hostname
date

# load software
module load stacks/2.53

# input, output files, directories

INDIR=../results/stacks/denovo

POPMAP=../meta/popmap_total.txt

populations \
-P $INDIR \
-M $POPMAP \
-p 1 \
-r 1 \
--hwe \
--genepop \
--vcf \
--fasta-samples \
--fasta-loci \
--treemix \
--structure \
-t 8