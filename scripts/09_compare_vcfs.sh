#!/bin/bash
#SBATCH --job-name=compare_vcfs
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

#####################################
# compare freebayes and stacks refmap
#####################################

module load bcftools/1.9
module load htslib/1.9
module load vt/0.57721
module load bedtools/2.29.0

INDIR=../results/filtered_vcfs

##############################
# check out filtered results:
##############################

vt peek $INDIR/fb_final.vcf.gz
vt peek $INDIR/denovo_final.vcf.gz
vt peek $INDIR/refmap_final.vcf.gz


##############################
# compare variant sets
##############################

# for two sets: VT partition
	# we can't directly compare the de novo set because it doesn't have genome coordinates
vt partition $INDIR/fb_final.vcf.gz $INDIR/refmap_final.vcf.gz

# extract intersections with bcftools isec
bcftools isec -O z -p $INDIR/isec $INDIR/fb_final.vcf.gz $INDIR/refmap_final.vcf.gz

#####################################
# can we investigate discrepancies?
#####################################

# sbf1.bed.gz gives the location of all sbf1 sites in the genome. 
# they need to be matched with mseI sites, so not all will be ddRAD sites
# we'll count up how many variants fall within 400bp of an sbf1 site
# we'll do this for stacks refmap and freebayes variant sets

# 'bedtools slop' creates a window of 400bp on either side of the sbf1 cut site
# 'bedtools map' counts the vcf records that fall in each sbf1 window

# bedtools wants to know the order of chromosomes/contigs in the vcf, this file provides that
FAI=../genome/GCA_023634145.1_fThyThy.pri.20220222_genomic.fna

bedtools slop -i ../meta/sbf1.bed.gz -l 400 -r 400 -g $FAI | \
bedtools map -a stdin -b ../results/filtered_vcfs/fb_final.vcf.gz -c 1 -o count -g $FAI \
>$INDIR/fb_sbf1_map.bed

bedtools slop -i ../meta/sbf1.bed.gz -l 400 -r 400 -g $FAI | \
bedtools map -a stdin -b ../results/filtered_vcfs/refmap_final.vcf.gz -c 1 -o count -g $FAI \
>$INDIR/refmap_sbf1_map.bed

# combine the counts files
	# col 4 = fb, col 4 = refmap
paste $INDIR/fb_sbf1_map.bed <(cut -f 4 $INDIR/refmap_sbf1_map.bed) >$INDIR/sbf1_vcf_counts.bed

# how many variants are sbf1-associated?
# almost all variants are associated with sbf1 sites in the reference genome. 
	# probably most others are in sites that have mutations in the ref genome, but not in most of the population

awk '{x+=$4}{y+=$5}END{print x " sbf1-associated variants in fb\n" y " sbf1-associated variants in stacks refmap\n"}' $INDIR/sbf1_vcf_counts.bed

# how many sbf1 sites have variants for each approach?
awk '{if ($4 > 0) x += 1}END{print x " sbf1 sites have variants in fb"}' $INDIR/fb_sbf1_map.bed

awk '{if ($4 > 0) x += 1}END{print x " sbf1 sites have variants in stacks refmap"}' $INDIR/refmap_sbf1_map.bed

# so about 6100 sbf1 sites don't have called variants in stacks refmap. why? we would have to explore more deeply. 

