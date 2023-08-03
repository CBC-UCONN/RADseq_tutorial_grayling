#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-552]%100

hostname
date

##################################
# align sequences to the reference
##################################

# this script aligns sequences to the reference genome

# load software--------------------------------------------------------------------------------
module load bwa/0.7.17
module load samtools/1.16.1

# input, output files and directories----------------------------------------------------------
INDIR=../data/demux

OUTDIR=../results/aligned
mkdir -p $OUTDIR

# reference genome index
REFERENCE=../genome/grayling

# make a bash array of fastq files
FASTQS=($(ls $INDIR/*1.fq.gz | grep -v "rem...fq.gz"))

# pull out a single fastq file pair
FQ1=$(echo ${FASTQS[$SLURM_ARRAY_TASK_ID]})
FQ2=$(echo $FQ1 | sed 's/1.fq.gz/2.fq.gz/')

# get sample ID
SAM=$(basename $FQ1 .1.fq.gz)

# create bam file name:
BAM=${SAM}.bam

# get the sample ID and use it to specify the read group. 
RG=$(echo \@RG\\tID:$SAM\\tSM:$SAM)

echo $OUTDIR
echo $BAM

# align sequences--------------------------------------------------------------------------------
# run bwa mem to align, then pipe it to samtools to compress, then again to sort
bwa mem -t 4 -R $RG $REFERENCE $FQ1 $FQ2 | \
samtools view -S -h -u - | \
samtools sort -T /scratch/${SAM}_${USER} - >$OUTDIR/$BAM

# index the bam file
samtools index $OUTDIR/$BAM

date

