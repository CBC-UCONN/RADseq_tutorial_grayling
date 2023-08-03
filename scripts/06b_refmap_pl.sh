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

##################################
# run stacks refmap.pl
##################################

# this script runs the stacks reference mapping module perl wrapper

# load software------------------------------------------------------------
module load stacks/2.64
module load GATK/4.1.8.1
module load bcftools/1.9
module load htslib/1.9

# input, output files, directories-----------------------------------------
INDIR=../results/aligned

OUTDIR=../results/stacks/refmap
mkdir -p $OUTDIR

# popmap file
POPMAP=../meta/popmap_total.txt

# run refmap.pl------------------------------------------------------------
	# provide `populations` options
ref_map.pl \
--samples $INDIR \
--popmap $POPMAP \
-o $OUTDIR \
-T 10 \
-X "populations:-p 1" \
-X "populations:-R 0.5" \
-X "populations:--genepop" \
-X "populations:--hwe" \
-X "populations:--vcf" \
-X "populations:--treemix" \
-X "populations:--structure" \
-X "populations:--fasta-samples" \
-X "populations:--fasta-loci"


# modify vcf file------------------------------------------------------------
# stacks doesn't include a sequence dictionary in the header.
# this is required by some tools. add it using gatk. 

ALIGNDIR=../results/aligned
BAMDICT=../results/aligned/Golden1A06.bam

gatk UpdateVCFSequenceDictionary \
     -V $OUTDIR/populations.snps.vcf \
     --source-dictionary $BAMDICT \
     --output $OUTDIR/populations.snps.dict.vcf \
     --replace true

# stacks does not properly sort the vcf file. sort it. 
bcftools sort --max-mem 30G -O z $OUTDIR/populations.snps.dict.vcf >$OUTDIR/populations.snps.dict.vcf.gz

# index the vcf file. 
tabix -p vcf $OUTDIR/populations.snps.dict.vcf.gz

