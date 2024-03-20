#!/bin/bash
#SBATCH --job-name=r80_optimization
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[1-64]%5

hostname
date

############################
# do r80 optimization
############################

# this is an array job. the code will be run 64 times, once for pair of parameters.
# we use the SLURM_ARRAY_TASK_ID to create the parameters file and pull out a single pair for each iteration 

# load software--------------------------------------------------------------------------------
module load stacks/2.64

#input/output directories, supplementary files-------------------------------------------------
INDIR=../data/demux

# make output directory if it doesn't exist
SUMMARYDIR=../results/stacks/denovo_r80/
mkdir -p ${SUMMARYDIR}

OUTDIR=../results/stacks/denovo_r80/iteration_${SLURM_ARRAY_TASK_ID}
mkdir -p ${OUTDIR}

# create file with parameter combinations by looping over all values of m and M
ALLm=($(echo {1..8}))
ALLM=($(echo {3..10}))

for m in ${ALLm[@]}
    do 
    for M in ${ALLM[@]}
        do 
        echo -e "${m}\t${M}"
    done
done >${OUTDIR}/params.txt

# select a pair of parameters
LITTLEM=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${OUTDIR}/params.txt | cut -f 1)
BIGM=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${OUTDIR}/params.txt | cut -f 2)

# use limited popmap to make this run faster
POPMAP=../meta/popmap_cstacks.txt

echo -e "PARAMETERS:\t${LITTLEM}\t${BIGM}"

denovo_map.pl \
    -T 12 \
    -o ${OUTDIR} \
    --popmap ${POPMAP} \
    --samples ${INDIR} \
    --paired \
    -X "ustacks:-m ${LITTLEM}" \
    -X "ustacks:-M ${BIGM}" \
    -X "cstacks:-n ${BIGM}" \
    -X "populations:-p 1" \
    -X "populations:-R 0.8" \
    -X "populations:--genepop" \
    -X "populations:--hwe" \
    -X "populations:--vcf" \
    -X "populations:--treemix" \
    -X "populations:--structure" \
    -X "populations:--fasta-samples" \
    -X "populations:--fasta-loci"

# to save space, we will simply grab the r80 statistic and then delete the entire run. It will have to be run again with selected parameters. 
    # we set -R = 0.8 to emit only loci present in 80% of individuals
    # then we simply count how many records are emitted in the haplotype vcf file (snps within a locus are concatenated into haplotypes in this file)
    # you don't need to re-run with -R 0.8. you can be more or less stringent depending on your aims. 
    
R80=$(grep -v "^#" ${OUTDIR}/populations.haps.vcf | wc -l)

echo -e "${LITTLEM}\t${BIGM}\t${R80}" >${SUMMARYDIR}/R80_${LITTLEM}_${BIGM}.txt

rm -r ${OUTDIR}
