#!/bin/bash
#SBATCH --job-name=ustacks
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-552]%20

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID
date

############################
# run `ustacks`
############################

# ustacks is the first step of the stacks de novo pipeline

# this is an array job. the code will be run 33 times, once for each sample, with 20 jobs running simultaneously.
# each instance will have a SLURM_ARRAY_TASK_ID between 0 and 32. 
# we use the SLURM_ARRAY_TASK_ID to grab a different sample for each job. 

# load software--------------------------------------------------------------------------------
module load stacks/2.64

#input/output directories, supplementary files-------------------------------------------------
INDIR=../data/demux

# make output directory if it doesn't exist
OUTDIR=../results/stacks/denovo
mkdir -p $OUTDIR

# make an array of all R1 fastq files, exclude rem.x.fq.gz, they're the failed reads
FASTQS=($(ls $INDIR/*1.fq.gz | grep -v "rem...fq.gz"))
# select the fastq file for this array task using SLURM_ARRAY_TASK_ID
INFILE=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]})

# specify the unique integer (-i $ID), add 1 b/c the task ids start at 0. 
ID=$(expr 1 + $SLURM_ARRAY_TASK_ID)
# specify the sample ID (--name $SAM)
# pull the sample ID from the fastq file name using grep
SAM=$(basename $INFILE .1.fq.gz)

# run ustacks--------------------------------------------------------------------------------
ustacks \
-f $INFILE \
-o $OUTDIR \
-i $ID \
--name $SAM \
-t gzfastq \
-p 6 \
-M 2 \
-m 3 \
