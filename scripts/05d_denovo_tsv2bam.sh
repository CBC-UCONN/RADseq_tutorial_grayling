#!/bin/bash 
#SBATCH --job-name=tsv2bam
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

# load software
module load stacks/2.53

# input, output files, directories

INDIR=../results/stacks/denovo
PEDIR=../data/demux
POPMAP=../meta/popmap_total.txt

tsv2bam \
-P $INDIR \
-M $POPMAP \
-t 10 \
-R $PEDIR

date