#!/bin/bash
#BSUB -J concat_impute
#BSUB -R "rusage[mem=182000]"
#BSUB -o concat_impute.out
#BSUB -e concat_impute.err
#BSUB -n 1
#BSUB -q short
#BSUB -W 02:00

module load bcftools
module load plink2/alpha6.1amd

# ---- PARAMETERS  ----
VCF_DIR=/pi/manuel.garber-umw/human/VIGOR/groberts/imputed_genotypes/combined_output
OUT_DIR=$VCF_DIR/final_output
mkdir -p "$OUT_DIR"

# Chromosome list (space-separated). Can include X, Y, MT if needed.
CHROM_LIST="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

# VCF filename prefix pattern (must include "chr" placeholder).
# E.g., for files like "chr1_imputed_combined_hg19.vcf.gz", use:
HG19_PREFIX="chr{CHR}_imputed_combined_hg19.vcf.gz"
HG38_PREFIX="chr{CHR}_imputed_combined_hg38.vcf.gz"

# Output file bases
HG19_OUT_BASE=VIGOR_imputed_hg19
HG38_OUT_BASE=VIGOR_imputed_hg38

# ---- Function: Expand file list for bcftools concat ----
build_vcf_list() {
  local pattern="$1"
  local list=""
  for CHR in $CHROM_LIST; do
    file="${VCF_DIR}/$(echo "$pattern" | sed "s/{CHR}/$CHR/")"
    if [[ -f "$file" ]]; then
      list+="$file "
    else
      echo "Warning: Missing file for chr$CHR: $file" >&2
    fi
  done
  echo "$list"
}

# ---- Step 1: Concatenate hg19 VCFs ----
echo "Concatenating hg19 VCFs..."
HG19_VCF_LIST=$(build_vcf_list "$HG19_PREFIX")
HG19_OUT=${OUT_DIR}/${HG19_OUT_BASE}.vcf.gz

bcftools concat -a -O z $HG19_VCF_LIST -o "$HG19_OUT"
bcftools index -f "$HG19_OUT"

# ---- Step 2: Concatenate hg38 VCFs ----
echo "Concatenating hg38 VCFs..."
HG38_VCF_LIST=$(build_vcf_list "$HG38_PREFIX")
HG38_OUT=${OUT_DIR}/${HG38_OUT_BASE}.vcf.gz

bcftools concat -a -O z $HG38_VCF_LIST -o "$HG38_OUT"
bcftools index -f "$HG38_OUT"
