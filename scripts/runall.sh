#!/bin/bash

jid1=$(sbatch --parsable 00_symlinkdata.sh )
jid2=$(sbatch --parsable --dependency=afterok:$jid1 01_get_genome.sh )
jid3=$(sbatch --parsable --dependency=afterok:$jid2 02_process_radtags.sh )
jid4=$(sbatch --parsable --dependency=afterok:$jid3 03_fastqc_raw.sh )
jid5=$(sbatch --parsable --dependency=afterok:$jid4 04_multiqc.sh )
jid6=$(sbatch --parsable --dependency=afterok:$jid5 05a_denovo_ustacks.sh )
jid7=$(sbatch --parsable --dependency=afterok:$jid6 05b_denovo_cstacks.sh )
jid8=$(sbatch --parsable --dependency=afterok:$jid7 05c_denovo_sstacks.sh )
jid9=$(sbatch --parsable --dependency=afterok:$jid8 05d_denovo_tsv2bam.sh )
jid10=$(sbatch --parsable --dependency=afterok:$jid9 05e_denovo_gstacks.sh )
jid11=$(sbatch --parsable --dependency=afterok:$jid10 05f_denovo_populations.sh )
jid12=$(sbatch --parsable --dependency=afterok:$jid11 06a_refmap_align.sh )
jid13=$(sbatch --parsable --dependency=afterok:$jid12 06b_refmap_pl.sh )
jid14=$(sbatch --parsable --dependency=afterok:$jid13 06c_alignQC.sh )
jid15=$(sbatch --parsable --dependency=afterok:$jid14 07_freebayes_parallel.sh )
jid16=$(sbatch --parsable --dependency=afterok:$jid15 08_filter_vcfs.sh )
jid17=$(sbatch --parsable --dependency=afterok:$jid16 09_compare_vcfs.sh )
jid18=$(sbatch --parsable --dependency=afterok:$jid17 10_reformat.sh)
jid19=$(sbatch --parsable --dependency=afterok:$jid18 11_analysis.sh)
jid20=$(sbatch --parsable --dependency=afterok:$jid19 05g_r80.sh)
