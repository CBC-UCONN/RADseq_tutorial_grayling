#!/bin/bash
#SBATCH --job-name=fastqc_raw
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --array=[0-552]
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

##########################
# run fastqc 
##########################

# load software
module load fastqc/0.11.7

#input/output directories, supplementary files
INDIR=../data/demux/

OUTDIR=../results/fastqc
mkdir -p $OUTDIR

# make bash array of R1 fastqs
FQ=($(find $INDIR -name "*1.fq.gz" | grep -v "rem...fq.gz" | sort))

# get R1 and R2 file names for this array instance
FQ1=${FQ[$SLURM_ARRAY_TASK_ID]}
FQ2=$(echo $FQ1 | sed 's/1.fq.gz/2.fq.gz/')

# run fastqc on read 1, read 2
fastqc -t 2 -o $OUTDIR $FQ1
fastqc -t 2 -o $OUTDIR $FQ2
