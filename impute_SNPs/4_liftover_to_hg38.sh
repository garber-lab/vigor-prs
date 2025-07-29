#!/bin/bash
#BSUB -J liftover_chr22
#BSUB -R "rusage[mem=128000]"
#BSUB -o liftover_chr22.out
#BSUB -e liftover_chr22.err
#BSUB -q short
#BSUB -W 0:30
#BSUB -n 1

module load picard/3.1.1
module load bcftools/1.16
module load htslib

# ---- PARAMETERS ----
CHAIN=/home/genevieve.roberts-umw/liftover_chain/hg19ToHg38.over.chain
REFERENCE=/home/genevieve.roberts-umw/liftover_chain/hg38.fa

IN_DIR=/home/genevieve.roberts-umw/imputed_genotypes/combined_output
OUT_DIR=/home/genevieve.roberts-umw/imputed_genotypes/combined_output
CHR=22
VCF_PREFIX=chr${CHR}_imputed_combined

# ---- FILE PATHS ----
IN_VCF=${IN_DIR}/${VCF_PREFIX}.vcf.gz
CLEANED_VCF=${OUT_DIR}/tmp_${VCF_PREFIX}_cleaned.vcf.gz
VCF_WITH_CHR=${OUT_DIR}/tmp_${VCF_PREFIX}_addchr.vcf.gz
REJECT_VCF=${OUT_DIR}/tmp_grch38_chr${CHR}_rejected.vcf.gz
LIFTED_VCF=${OUT_DIR}/grch38_chr${CHR}.vcf.gz

# ---- STEP 1: Remove symbolic reference alleles ----
echo "Filtering out symbolic reference alleles like <CN0> from REF field..."
bcftools view -e 'REF ~ "^<"' "$IN_VCF" -O z -o "$CLEANED_VCF"
bcftools index -f "$CLEANED_VCF"

# ---- STEP 2: Add 'chr' prefix to contig ----
echo "Adding 'chr' prefix to contig ${CHR}..."
bcftools annotate \
  --rename-chrs <(echo -e "${CHR}\tchr${CHR}") \
  -O z -o "$VCF_WITH_CHR" "$CLEANED_VCF"
bcftools index -f "$VCF_WITH_CHR"

# ---- STEP 3: Run Picard LiftoverVcf ----
echo "Running Picard LiftoverVcf..."
JAVA_MEM=120g
PICARD_OPTIONS="-Xmx$JAVA_MEM"

java $PICARD_OPTIONS -jar $PICARDJAR LiftoverVcf \
  I="$VCF_WITH_CHR" \
  O="$LIFTED_VCF" \
  CHAIN="$CHAIN" \
  REJECT="$REJECT_VCF" \
  R="$REFERENCE" \
  MAX_RECORDS_IN_RAM=1000000 \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT

# Remove temporary files
rm -f "$CLEANED_VCF" "$VCF_WITH_CHR" "$REJECT_VCF"

echo "Liftover complete."
