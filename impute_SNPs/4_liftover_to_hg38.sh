#!/bin/bash
#BSUB -J liftover[1-22]
#BSUB -R "rusage[mem=182000]"
#BSUB -o liftover_chr%I.out
#BSUB -e liftover_chr%I.err
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1

module load picard/3.1.1
module load bcftools/1.16
module load htslib

# ---- PARAMETERS ----
CHAIN=/home/genevieve.roberts-umw/liftover_chain/hg19ToHg38.over.chain
REFERENCE=/home/genevieve.roberts-umw/liftover_chain/hg38.fa

IN_DIR=/home/genevieve.roberts-umw/imputed_genotypes/combined_output
OUT_DIR=/home/genevieve.roberts-umw/imputed_genotypes/combined_output

CHR=${LSB_JOBINDEX}
VCF_PREFIX=chr${CHR}_imputed_combined_hg19

# ---- FILE PATHS ----
IN_VCF=${IN_DIR}/${VCF_PREFIX}.vcf.gz
REJECT_VCF=${OUT_DIR}/chr${CHR}_liftover_rejected.vcf.gz
LIFTED_VCF=${OUT_DIR}/chr${CHR}_imputed_combined_hg38.vcf.gz
ID_MAP_TXT_RAW=${OUT_DIR}/chr${CHR}_id_mapping.txt

# ---- STEP 1: Run Picard LiftoverVcf ----
echo "Running Picard LiftoverVcf for chr${CHR}..."
JAVA_MEM=180g
PICARD_OPTIONS="-Xmx$JAVA_MEM"

java $PICARD_OPTIONS -jar $PICARDJAR LiftoverVcf \
  I="$IN_VCF" \
  O="$LIFTED_VCF" \
  CHAIN="$CHAIN" \
  REJECT="$REJECT_VCF" \
  R="$REFERENCE" \
  MAX_RECORDS_IN_RAM=1000000 \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  RECOVER_SWAPPED_REF_ALT=true

rm -f "$REJECT_VCF"

# ---- STEP 2: Clean header (remove ##contig lines) ----
bcftools view -h "$LIFTED_VCF" | grep -v "^##contig=" > tmp_header.txt
bcftools reheader -h tmp_header.txt -o "${LIFTED_VCF}.tmp" "$LIFTED_VCF"
mv "${LIFTED_VCF}.tmp" "$LIFTED_VCF"
rm tmp_header.txt

# ---- STEP 3: Rename variant IDs based on new position ----
echo "Renaming variant IDs for chr${CHR}..."

# Extract header
bcftools view -h "$LIFTED_VCF" > vcf_header.tmp

# Rewrite body with updated IDs and capture ID mapping
bcftools view -H "$LIFTED_VCF" | \
  awk -v chr="$CHR" -v mapfile="$ID_MAP_TXT_RAW" 'BEGIN{OFS="\t"} {
    old_id = $3;
    if (old_id ~ /^HLA/) {
        new_id = old_id;
    } else {
        new_id = chr ":" $2 ":" $4 ":" $5;
    }
    $3 = new_id;
    print old_id, new_id >> mapfile;
    print $0;
  }' > body.tmp.vcf

# Compress the ID mapping file
gzip "$ID_MAP_TXT_RAW"

# Combine header and modified body, compress and index
cat vcf_header.tmp body.tmp.vcf | bgzip > "$LIFTED_VCF"
tabix -p vcf "$LIFTED_VCF"

# Clean up
rm vcf_header.tmp body.tmp.vcf
