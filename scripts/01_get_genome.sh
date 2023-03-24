#!/bin/bash
#SBATCH --job-name=get_genome
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

######################
# get reference genome
######################

# download genome assembly for grayling, Thymallus thymallus

# output directory
GENOMEDIR=../genome
mkdir -p $GENOMEDIR

# download genome
wget \
-P $GENOMEDIR \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/634/145/GCA_023634145.1_fThyThy.pri.20220222/GCA_023634145.1_fThyThy.pri.20220222_genomic.fna.gz

# decompress the genome
gunzip $GENOMEDIR/GCA_023634145.1_fThyThy.pri.20220222_genomic.fna.gz


# index the genome using bwa
module load bwa/0.7.17
bwa index \
-p $GENOMEDIR/grayling \
$GENOMEDIR/GCA_023634145.1_fThyThy.pri.20220222_genomic.fna

# index the genome using samtools
module load samtools/1.10
samtools faidx $GENOMEDIR/GCA_023634145.1_fThyThy.pri.20220222_genomic.fna


# run an R script to find all sbf1 sites, put them in a bed file in ../meta	
	# precomputed. output found in directory ../meta

# module load R/3.6.3

# Rscript find_sbf1_sites.R
# gzip ../meta/sbf1.bed
# gzip ../meta/sbf1_off.bed
