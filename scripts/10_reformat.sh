#!/bin/bash
#SBATCH --job-name=reformat
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=xeon
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#####################################
# reformat data
#####################################

module load bcftools/1.9
module load plink/2.00a2.3LM

# note that plink2 requires avx instructions (I think?) which requires a xeon processor
	# submit to xeon partition

# input/output directories
INDIR=../results/filtered_vcfs
OUTDIR=../results/variants_reformatted
mkdir -p $OUTDIR

# generate plink formatted file
plink2 --vcf $INDIR/fb_final.vcf.gz --allow-extra-chr --make-bed --out $OUTDIR/fb 
	# change chromosome names to integers so ADMIXTURE doesn't complain
	sed -i -E 's/CM(......).1/\1/; s/QMII(........).1/\1/' $OUTDIR/fb.bim

# create a reformatted file with genotypes only
bcftools query -H -f '%CHROM\t%POS\t[\t%GT]\n' $INDIR/fb_final.vcf.gz >$OUTDIR/fb_genos.tsv

bcftools query -H -f '%CHROM\t%POS\t[\t%GT]\n' $INDIR/refmap_final.vcf.gz >$OUTDIR/refmap_genos.tsv



