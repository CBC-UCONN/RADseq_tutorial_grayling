#!/bin/bash
#SBATCH --job-name=raw_data_symlinks
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=128M
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err


#####################
# symlink data
#####################

# this script symlinks the raw data pools to the working directory
	# instead of copying the raw data, we create pointers to it. 

# input, output directories
POOLPATH="/core/labs/Wegrzyn/Urban_RAD_ArcticGrayling/raw_fastq/raw_fastq_files"

RAWDATADIR=../data/pools
mkdir -p $RAWDATADIR

# using a for loop, create a symlink for each fastq.gz file
for f in ${POOLPATH}/Golden-Pool*fastq.gz; do
        # write filename and file base name to standard out
        echo $f
        echo $(basename ${f})

        # symlink file	
        ln -s ${f} $RAWDATADIR/$(basename ${f})
done

