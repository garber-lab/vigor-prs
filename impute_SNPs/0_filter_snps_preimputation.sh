#!/bin/bash
#BSUB -J pre_imputation_filt
#BSUB -R "rusage[mem=4096]"
#BSUB -o pre_imputation_filt.out
#BSUB -e pre_imputation_filt.err
#BSUB -q short
#BSUB -n 1

# Load necessary modules
module load plink2
module load htslib

# Ensure output directory exists
mkdir -p /home/genevieve.roberts-umw/pre_imputation_vcf

for chr in {1..23}; do
  plink \
    --bfile "/pi/manuel.garber-umw/human/VIGOR/tj/ibd/run_1/all_samples_ibd" \
    --chr $chr \
    --hwe .001 \
    --maf .01 \
    --recode vcf \
    --out "/home/genevieve.roberts-umw/pre_imputation_vcf/chr${chr}_all_samples_ibd_hwe_0.001"

  VCF_FILE="/home/genevieve.roberts-umw/pre_imputation_vcf/chr${chr}_all_samples_ibd_hwe_0.001.vcf"

  if [[ -f "$VCF_FILE" ]]; then
    # Update chromosome codes
    awk '
      BEGIN { OFS="\t" }
      !/^##/ {
        if ($1 ~ /^[0-9]+$/) $1 = "chr" $1;
        else if ($1 == "23") $1 = "chrX";
        else if ($1 == "24") $1 = "chrY";
        else if ($1 == "26") $1 = "chrMT";
        print
      }
      /^##/ { print }
    ' "$VCF_FILE" > "${VCF_FILE}.tmp" && mv "${VCF_FILE}.tmp" "$VCF_FILE"

    # bgzip the updated VCF
    bgzip "$VCF_FILE"
  else
    echo "Error: Missing output file for chr${chr}" >&2
  fi
done

# Clean up logs
mkdir -p /home/genevieve.roberts-umw/pre_imputation_vcf/logs
mv /home/genevieve.roberts-umw/pre_imputation_vcf/*.log /home/genevieve.roberts-umw/pre_imputation_vcf/logs 2>/dev/null
rm -f /home/genevieve.roberts-umw/pre_imputation_vcf/*.hh
rm -f /home/genevieve.roberts-umw/pre_imputation_vcf/*.nosex
