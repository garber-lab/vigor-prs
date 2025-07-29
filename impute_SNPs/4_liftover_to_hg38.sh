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

# ---- PARAMETERS ----
CHAIN=/home/genevieve.roberts-umw/liftover_chain/hg19ToHg38.over.chain
REFERENCE=/home/genevieve.roberts-umw/liftover_chain/hg38.fa

IN_DIR=/home/genevieve.roberts-umw/imputed_genotypes/1KG
OUT_DIR=/home/genevieve.roberts-umw/imputed_genotypes/1KG
CHR=22
VCF_PREFIX=test_chr${CHR}

# ---- FILE PATHS ----
IN_VCF=${IN_DIR}/${VCF_PREFIX}.dose.vcf
VCF_1000=${OUT_DIR}/${VCF_PREFIX}_first1000.dose.vcf
OUT_VCF_WITH_CHR=${OUT_DIR}/${VCF_PREFIX}_first1000_addchr.dose.vcf.gz
REJECT_VCF=${OUT_DIR}/test_hg38_chr${CHR}_first1000_rejected.vcf
LIFTED_VCF=${OUT_DIR}/test_hg38_chr${CHR}_first1000.dose.vcf

# ---- STEP 0: Subset first 1000 SNPs ----
echo "Subsetting first 1000 records for chr${CHR}..."
(bcftools view -h "$IN_VCF" && bcftools view -H "$IN_VCF" | head -n 1000) > "$VCF_1000"

# ---- STEP 1: Add chr prefixes to contigs ----
echo "Adding 'chr' prefixes to contigs..."
bcftools annotate \
  --rename-chrs <(echo -e "${CHR}\tchr${CHR}") \
  -O z -o "$OUT_VCF_WITH_CHR" "$VCF_1000"

bcftools index "$OUT_VCF_WITH_CHR"

# ---- STEP 2: Run Picard LiftoverVcf ----
echo "Running LiftoverVcf..."
JAVA_MEM=120g
java -Xmx$JAVA_MEM -jar $PICARDJAR LiftoverVcf \
    I="$OUT_VCF_WITH_CHR" \
    O="$LIFTED_VCF" \
    CHAIN="$CHAIN" \
    REJECT="$REJECT_VCF" \
    R="$REFERENCE" \
    MAX_RECORDS_IN_RAM=1000000

echo "Liftover complete."
