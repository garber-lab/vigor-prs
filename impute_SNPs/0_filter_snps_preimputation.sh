#!/bin/bash
#BSUB -J pre_imputation_filt
#BSUB -R rusage[mem=4096]
#BSUB -o pre_imputation_filt.out
#BSUB -e pre_imputation_filt.err
#BSUB -q short
#BSUB -n 1

# Load necessary modules
module load plink2

#filter based on HWE
for chr in {1..23}; do
  plink \
    --bfile "/pi/manuel.garber-umw/human/VIGOR/tj/ibd/run_1/all_samples_ibd" \
    --chr $chr \
    --hwe .001 \
    --maf .01 \
    --recode vcf \
    --out "/home/genevieve.roberts-umw/pre_imputation_vcf/chr${chr}_all_samples_ibd_hwe_0.001"
done
