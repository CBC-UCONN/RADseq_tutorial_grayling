#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load PopLDdecay/3.41
module load plink/2.00a2.3LM
module load admixture/1.3.0

#######################################
# input/output directories
#######################################

SCRIPTDIR=$(pwd)

INDIR=../results/variants_reformatted

ADMIXDIR=../results/admixture
mkdir -p $ADMIXDIR/fb

PCADIR=../results/plink_pca
mkdir -p $PCADIR

#######################################
# population structure using admixture
#######################################

cd $ADMIXDIR
for K in {1..10}; \
	do admixture --cv $INDIR/fb.bed $K | tee log${K}.out
done
cd $SCRIPTDIR

#######################################
# pca using plink
#######################################

plink2 --bfile $INDIR/fb --pca --out $PCADIR/fb --allow-extra-chr

#######################################

#######################################


# fst using plink

	# needs the newest version...
