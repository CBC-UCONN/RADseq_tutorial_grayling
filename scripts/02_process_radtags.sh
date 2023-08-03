#!/bin/bash
#SBATCH --job-name=process_radtags
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=40G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --array=[0-11]
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err


############################
# demultiplex pooled data
############################

# This script demultiplexes each pool of RAD-seq data
# It's an array script, so it submits 12 jobs, one for each pool

module load stacks/2.64

# set input/output directories
INDIR=../data/pools/
META=../meta

OUTDIR=../data/demux/
mkdir -p $OUTDIR

# make a bash array of Read 1 pool fastq files
POOL=($(ls $INDIR/*R1*))

# make a bash array of 3-digit indexes for each pool:
IND=({001..012})

# get fastq and barcode file names for this array instance:
	# variable $SLURM_ARRAY_TASK_ID changes for each array instance

FQ1=${POOL[$SLURM_ARRAY_TASK_ID]}
FQ2=$(echo ${POOL[$SLURM_ARRAY_TASK_ID]} | sed 's/R1/R2/')

IND1=${IND[$SLURM_ARRAY_TASK_ID]}
BARCODE=$META/barcodes_${IND1}.txt

# run stacks process_radtags
process_radtags \
-i gzfastq \
-1 $FQ1 \
-2 $FQ2 \
-b $BARCODE \
-o $OUTDIR \
-c -q -r \
--inline_null \
--renz_1 sbfI --renz_2 mseI
