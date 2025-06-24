#!/bin/bash
#BSUB -J pre_imputation_filt
#BSUB -R "rusage[mem=4096]"
#BSUB -o pre_imputation_filt.out
#BSUB -e pre_imputation_filt.err
#BSUB -q short
#BSUB -n 1

# Load necessary modules
module load plink2/alpha6.1amd
module load htslib

# Ensure output directory exists
mkdir -p /home/genevieve.roberts-umw/pre_imputation_vcf

# Define paths to SNP list files
STARTING_BFILE="/pi/manuel.garber-umw/human/VIGOR/tj/ibd/run_2/all_samples_ibd"
SNPS_TO_REMOVE="/home/genevieve.roberts-umw/vigor-prs/impute_SNPs/drop_invalid_allele_snps.txt"
SNPS_TO_SET_REF="/home/genevieve.roberts-umw/vigor-prs/impute_SNPs/combined_ref_alleles.txt"
MHC_SNPS_TO_SET_REF="/home/genevieve.roberts-umw/vigor-prs/impute_SNPs/MHC_set_ref_snps.txt"
INDIVS_TO_REMOVE="/home/genevieve.roberts-umw/vigor-prs/impute_SNPs/sex_problems.txt"

for chr in {1..23}; do
  plink2 \
    --bfile "$STARTING_BFILE" \
    --chr $chr \
    --hwe .001 \
    --maf .01 \
    --exclude "$SNPS_TO_REMOVE" \
    --remove "$INDIVS_TO_REMOVE" \
    --snps-only just-acgt \
    --rm-dup exclude-all \
    --ref-allele "$SNPS_TO_SET_REF" \
    --recode vcf \
    --out "/home/genevieve.roberts-umw/pre_imputation_vcf/chr${chr}_all_samples_ibd_hwe_0.001"

  VCF_FILE="/home/genevieve.roberts-umw/pre_imputation_vcf/chr${chr}_all_samples_ibd_hwe_0.001.vcf"

  if [[ -f "$VCF_FILE" ]]; then
    # Update chromosome codes
    awk '
      BEGIN { OFS="\t" }
      /^##fileformat=VCFv4.3/ { print "##fileformat=VCFv4.2"; next }
      !/^##/ {
        if ($1 ~ /^[0-9]+$/) $1 = "chr" $1;
        else if ($1 == "chr23" || $1 == "X" || $1 == "23") $1 = "chrX";
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

# Additional processing for MHC region using MHC-specific SNP ref file
plink2 \
  --bfile "$STARTING_BFILE" \
  --chr 6 \
  --hwe .001 \
  --maf .01 \
  --exclude "$SNPS_TO_REMOVE" \
  --snps-only just-acgt \
  --rm-dup exclude-all \
  --ref-allele "$MHC_SNPS_TO_SET_REF" \
  --recode vcf \
  --out "/home/genevieve.roberts-umw/pre_imputation_vcf/chr6_MHC_all_samples_ibd_hwe_0.001"

MHC_VCF="/home/genevieve.roberts-umw/pre_imputation_vcf/chr6_MHC_all_samples_ibd_hwe_0.001.vcf"

if [[ -f "$MHC_VCF" ]]; then
  awk '
      BEGIN { OFS="\t" }
      /^##fileformat=VCFv4.3/ { print "##fileformat=VCFv4.2"; next }
      !/^##/ {
      if ($1 ~ /^[0-9]+$/) $1 = "chr" $1;
      else if ($1 == "chr23" || $1 == "X" || $1 == "23") $1 = "chrX";
      else if ($1 == "24") $1 = "chrY";
      else if ($1 == "26") $1 = "chrMT";
      print
    }
    /^##/ { print }
  ' "$MHC_VCF" > "${MHC_VCF}.tmp" && mv "${MHC_VCF}.tmp" "$MHC_VCF"

  bgzip "$MHC_VCF"
else
  echo "Error: Missing MHC output file for chr6" >&2
fi

# Clean up logs
mkdir -p /home/genevieve.roberts-umw/pre_imputation_vcf/logs
mv /home/genevieve.roberts-umw/pre_imputation_vcf/*.log /home/genevieve.roberts-umw/pre_imputation_vcf/logs 2>/dev/null
rm -f /home/genevieve.roberts-umw/pre_imputation_vcf/*.hh
rm -f /home/genevieve.roberts-umw/pre_imputation_vcf/*.nosex
