#!/bin/bash
#SBATCH --job-name=filter_vcfs
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

hostname
date

##################################
# filter variant sets
##################################

module load bcftools/1.9
module load htslib/1.9
module load vcftools/0.1.16
module load vcflib/1.0.0-rc1


###############################
# set input, output directories
###############################

OUTDIR=../results/filtered_vcfs
mkdir -p $OUTDIR

DENOVO=../results/stacks/denovo
REFMAP=../results/stacks/refmap
FREEBA=../results/freebayes

#############################
# filter SITES by missingness
#############################

# also remove multiallelic sites and indels

# stacks de novo-----------------------------------------------------------------------
vcftools --gzvcf $DENOVO/populations.snps.vcf \
	--max-missing-count 116 --mac 3 --remove-indels --max-alleles 2 --min-alleles 2 \
	--recode \
	--stdout | \
	bgzip >$OUTDIR/stacks_denovo.vcf.gz

	# output missing individual report
	vcftools --gzvcf $OUTDIR/stacks_denovo.vcf.gz --out $OUTDIR/stacks_denovo --missing-indv

# stacks refmap-----------------------------------------------------------------------
vcftools --gzvcf $REFMAP/populations.snps.dict.vcf.gz \
	--max-missing-count 116 --mac 3 --remove-indels --max-alleles 2 --min-alleles 2 \
	--recode \
	--stdout | \
	bgzip >$OUTDIR/stacks_refmap.vcf.gz

	# output missing individual report
	vcftools --gzvcf $OUTDIR/stacks_refmap.vcf.gz --out $OUTDIR/stacks_refmap --missing-indv

# freebayes----------------------------------------------------------------------------
	# - note also that:
		# 1. bcftools norm normalizes variant representation
		# 2. vcfallelic primitives breaks down haplotype alleles into constituent parts
GEN=../genome/GCA_023634145.1_fThyThy.pri.20220222_genomic.fna
bcftools norm -f $GEN $FREEBA/fb_parallel.vcf.gz | \
	vcfallelicprimitives --keep-info --keep-geno | \
	vcfstreamsort | \
	vcftools --vcf - \
	--max-missing-count 116 --mac 3 --remove-indels --max-alleles 2 --min-alleles 2 \
	--recode \
	--stdout | bgzip >$OUTDIR/fb.vcf.gz

	# output missing individual report
	vcftools --gzvcf $OUTDIR/fb.vcf.gz --out $OUTDIR/fb --missing-indv

###################################
# filter INDIVIDUALS by missingness
###################################

# create a list of samples with high rates of missing genotypes to exclude
DROPSAMPLES=$(tail -n +2 $OUTDIR/fb.imiss | awk '$5 > .1' | cut -f 1 | tr "\n" "," | sed 's/,$//')

# drop samples with high rates of missing data, also exclude low variant quality variants from freebayes

# stacks denovo
bcftools view -s ^$DROPSAMPLES $OUTDIR/stacks_denovo.vcf.gz | bgzip -c > $OUTDIR/denovo_final.vcf.gz

# stacks refmap
bcftools view -s ^$DROPSAMPLES $OUTDIR/stacks_refmap.vcf.gz | bgzip -c > $OUTDIR/refmap_final.vcf.gz

# freebayes
bcftools view -s ^$DROPSAMPLES $OUTDIR/fb.vcf.gz | bcftools filter -s LowQual -e '%QUAL<50' | bgzip -c > $OUTDIR/fb_final.vcf.gz

##############
# make indexes
##############

for file in $OUTDIR/*vcf.gz
do tabix -f -p vcf $file
done
