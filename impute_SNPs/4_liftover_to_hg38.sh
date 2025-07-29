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
VCF_1000=${OUT_DIR}/${VCF_PREFIX}_first1000.vcf
VCF_1000_GZ=${VCF_1000}.gz
OUT_VCF_WITH_CHR=${OUT_DIR}/${VCF_PREFIX}_first1000_addchr.dose.vcf.gz
REJECT_VCF=${OUT_DIR}/test_hg38_chr${CHR}_first1000_rejected.vcf.gz
LIFTED_VCF=${OUT_DIR}/test_hg38_chr${CHR}_first1000.dose.vcf.gz

# ---- STEP 0: Subset first 1000 records for chr${CHR} ----
echo "Subsetting first 1000 records for chr${CHR}..."
# Save as uncompressed first
(bcftools view -h "$IN_VCF" && bcftools view -H "$IN_VCF" | head -n 1000) > "$VCF_1000"
bgzip -f "$VCF_1000"  # compress to .vcf.gz
tabix -f "$VCF_1000_GZ"  # index

# ---- STEP 1: Add 'chr' prefix to contig name ----
echo "Adding 'chr' prefix to contig ${CHR}..."
bcftools annotate \
  --rename-chrs <(echo -e "${CHR}\tchr${CHR}") \
  -O z -o "$OUT_VCF_WITH_CHR" "$VCF_1000_GZ"
bcftools index -f "$OUT_VCF_WITH_CHR"

# ---- STEP 2: Run Picard LiftoverVcf ----
echo "Running Picard LiftoverVcf..."
JAVA_MEM=120g
PICARD_OPTIONS="-Xmx$JAVA_MEM"

# Run liftover
java $PICARD_OPTIONS -jar $PICARDJAR LiftoverVcf \
  I="$OUT_VCF_WITH_CHR" \
  O="$LIFTED_VCF" \
  CHAIN="$CHAIN" \
  REJECT="$REJECT_VCF" \
  R="$REFERENCE" \
  MAX_RECORDS_IN_RAM=1000000 \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT

echo "Liftover complete."
