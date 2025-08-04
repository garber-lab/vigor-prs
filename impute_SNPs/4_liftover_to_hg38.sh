#!/bin/bash
#BSUB -J liftover[6]
#BSUB -R "rusage[mem=182000]"
#BSUB -o liftover_chr%I.out
#BSUB -e liftover_chr%I.err
#BSUB -q short
#BSUB -W 0:30
#BSUB -n 1

module load picard/3.1.1
module load bcftools/1.16
module load htslib

# ---- PARAMETERS ----
CHAIN=/pi/manuel.garber-umw/human/VIGOR/groberts/liftover_chain/hg19ToHg38.over.chain
REFERENCE=/pi/manuel.garber-umw/human/VIGOR/groberts/liftover_chain/hg38.fa
HG19_HLA_BED=/pi/manuel.garber-umw/human/VIGOR/groberts/liftover_chain/hla_hg19.bed
HG38_HLA_BED=/pi/manuel.garber-umw/human/VIGOR/groberts/liftover_chain/hla_hg38.bed

IN_DIR=/pi/manuel.garber-umw/human/VIGOR/groberts/imputed_genotypes/combined_output
OUT_DIR=/pi/manuel.garber-umw/human/VIGOR/groberts/imputed_genotypes/combined_output

CHR=${LSB_JOBINDEX}
VCF_PREFIX=chr${CHR}_imputed_combined_hg19

IN_VCF=${IN_DIR}/${VCF_PREFIX}.vcf.gz
REJECT_VCF=${OUT_DIR}/chr${CHR}_liftover_rejected.vcf.gz
LIFTED_VCF=${OUT_DIR}/chr${CHR}_imputed_combined_hg38.vcf.gz
ID_MAP_TXT_RAW=${OUT_DIR}/chr${CHR}_id_mapping.txt

# ---- STEP 1: Liftover using Picard ----
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

# ---- STEP 2: Clean header ----
bcftools view -h "$LIFTED_VCF" | grep -v "^##contig=" > "${OUT_DIR}/tmp_header_chr${CHR}.txt"
bcftools reheader -h "${OUT_DIR}/tmp_header_chr${CHR}.txt" -o "${LIFTED_VCF}.tmp" "$LIFTED_VCF"
mv "${LIFTED_VCF}.tmp" "$LIFTED_VCF"
rm -f "${OUT_DIR}/tmp_header_chr${CHR}.txt"

# ---- STEP 3: Rename variant IDs ----
echo "Renaming variant IDs for chr${CHR}..."
bcftools view -h "$LIFTED_VCF" > "${OUT_DIR}/vcf_header_chr${CHR}.tmp"
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
  }' > "${OUT_DIR}/body_chr${CHR}.tmp.vcf"

gzip "$ID_MAP_TXT_RAW"
cat "${OUT_DIR}/vcf_header_chr${CHR}.tmp" "${OUT_DIR}/body_chr${CHR}.tmp.vcf" | bgzip > "$LIFTED_VCF"
tabix -f -p vcf "$LIFTED_VCF"

rm -f "${OUT_DIR}/vcf_header_chr${CHR}.tmp" "${OUT_DIR}/body_chr${CHR}.tmp.vcf"

# ---- STEP 4: Handle HLA alleles for chr6 only ----
if [[ "$CHR" -eq 6 ]]; then
  echo "Handling HLA alleles for chr6..."

  # 1. Remove HLA variants from lifted VCF
  bcftools view -e 'ID ~ "^HLA_"' "$LIFTED_VCF" -O z -o "${OUT_DIR}/chr${CHR}_no_hla.vcf.gz"
  bcftools index -f "${OUT_DIR}/chr${CHR}_no_hla.vcf.gz"

  # 2. Extract HLA variants from original VCF
  bcftools view -R "$HG19_HLA_BED" "$IN_VCF" -O z -o "${OUT_DIR}/chr${CHR}_hla.vcf.gz"
  bcftools index -f "${OUT_DIR}/chr${CHR}_hla.vcf.gz"

  # 3. Save header separately
  bcftools view -h "${OUT_DIR}/chr${CHR}_hla.vcf.gz" > "${OUT_DIR}/hla_header_chr${CHR}.tmp.vcf"

  # 4. Update positions using BED liftover info
  awk 'BEGIN{OFS="\t"} FNR==NR {pos[$4]=$2; chr[$4]=$1; next}
       {
         id = $3;
         if (id in pos) {
           $1 = chr[id];
           $2 = pos[id];
         }
         print;
       }' "$HG38_HLA_BED" <(bcftools view -H "${OUT_DIR}/chr${CHR}_hla.vcf.gz") > "${OUT_DIR}/hla_lifted_body_chr${CHR}.vcf"

  # 5. Combine header + updated body
  cat "${OUT_DIR}/hla_header_chr${CHR}.tmp.vcf" "${OUT_DIR}/hla_lifted_body_chr${CHR}.vcf" | \
  bcftools sort -O z -o "${OUT_DIR}/chr${CHR}_hla_lifted.vcf.gz"
  bcftools index -f "${OUT_DIR}/chr${CHR}_hla_lifted.vcf.gz"
  rm -f "${OUT_DIR}/hla_header_chr${CHR}.tmp.vcf" "${OUT_DIR}/hla_lifted_body_chr${CHR}.vcf"

  # 6. Ensure sample names match
  bcftools reheader -s <(bcftools query -l "${OUT_DIR}/chr${CHR}_no_hla.vcf.gz") \
    -o "${OUT_DIR}/fixed_hla_chr${CHR}.vcf.gz" "${OUT_DIR}/chr${CHR}_hla_lifted.vcf.gz"
  bcftools index -f "${OUT_DIR}/fixed_hla_chr${CHR}.vcf.gz"

  # 7. Concatenate non-HLA + HLA
  bcftools concat -a -O z \
    -o "${OUT_DIR}/chr${CHR}_with_hla.vcf.gz" \
    "${OUT_DIR}/chr${CHR}_no_hla.vcf.gz" \
    "${OUT_DIR}/fixed_hla_chr${CHR}.vcf.gz"
  bcftools index -f "${OUT_DIR}/chr${CHR}_with_hla.vcf.gz"

  # Replace original lifted VCF
  mv "${OUT_DIR}/chr${CHR}_with_hla.vcf.gz" "$LIFTED_VCF"
  bcftools index -f "$LIFTED_VCF"

  # redo the ID mapping to include HLA variants
  rm -f "$ID_MAP_TXT_RAW.gz"
  bcftools view -h "$LIFTED_VCF" > "${OUT_DIR}/vcf_header_chr${CHR}.tmp"
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
    }' > "${OUT_DIR}/body_chr${CHR}.tmp.vcf"
  gzip "$ID_MAP_TXT_RAW"

  # 8. Cleanup
  rm -f \
    "${OUT_DIR}/chr${CHR}_no_hla.vcf.gz" \
    "${OUT_DIR}/chr${CHR}_no_hla.vcf.gz.tbi" \
    "${OUT_DIR}/chr${CHR}_hla.vcf.gz" \
    "${OUT_DIR}/chr${CHR}_hla.vcf.gz.tbi" \
    "${OUT_DIR}/chr${CHR}_hla_lifted.vcf.gz" \
    "${OUT_DIR}/chr${CHR}_hla_lifted.vcf.gz.tbi" \
    "${OUT_DIR}/chr${CHR}_with_hla.vcf.gz.tbi" \
    "${OUT_DIR}/fixed_hla_chr${CHR}.vcf.gz" \
    "${OUT_DIR}/fixed_hla_chr${CHR}.vcf.gz.tbi" \
    "${OUT_DIR}/vcf_header_chr${CHR}.tmp" \
    "${OUT_DIR}/body_chr${CHR}.tmp.vcf"
fi
