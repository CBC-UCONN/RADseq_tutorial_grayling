#!/bin/bash 
#SBATCH --job-name=refmap.pl
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --qos=general
#SBATCH --partition=general

hostname
date

# load software
module load stacks/2.53

# input, output files, directories
INDIR=../results/aligned

OUTDIR=../results/stacks/refmap
mkdir -p $OUTDIR

# popmap file
POPMAP=../meta/popmap_total.txt

# refmap.pl -s option is broken. 
ref_map.pl \
--samples $INDIR \
--popmap $POPMAP \
-o $OUTDIR \
-T 10 \
-X "populations:-p 1" \
-X "populations:-r 1" \
-X "populations:--genepop" \
-X "populations:--hwe" \
-X "populations:--vcf" \
-X "populations:--treemix" \
-X "populations:--structure" \
-X "populations:--fasta-samples" \
-X "populations:--fasta-loci"