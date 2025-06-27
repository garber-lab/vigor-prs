#!/bin/bash
#BSUB -J combine_imp_vcfs
#BSUB -R "rusage[mem=8192]"
#BSUB -o combine_imp_vcfs.out
#BSUB -e combine_imp_vcfs.err
#BSUB -q short
#BSUB -n 1

# Load necessary modules
module load htslib

#change working directory
cd /home/genevieve.roberts-umw/imputed_genotypes

# Output directory
OUTPUT_DIR="./combined_output"
mkdir -p "$OUTPUT_DIR"

#make a temporary directory
TMP_DIR="/home/genevieve.roberts-umw/imputed_genotypes/imputation_tmp"
mkdir -p "$TMP_DIR"

# Imputed sources
INPUT_DIRS=("1KG" "HRC" "MHC_alleles")

# Loop over chromosomes
for chr in {21..22}; do
    echo "Processing chromosome $chr..."

    ANNOTATED_VCFS=()

    # Step 1: Create site-only VCFs from each panel with SOURCE tag
    for src in "${INPUT_DIRS[@]}"; do
        INPUT_VCF="${src}/chr${chr}.dose.vcf.gz"
        if [[ -f "$INPUT_VCF" ]]; then
            SITE_ONLY="${TMP_DIR}/${src}_chr${chr}_site.vcf.gz"
            MODIFIED_SITE="${TMP_DIR}/${src}_chr${chr}_site_modified.vcf.gz"
            TAGGED_VCF="${TMP_DIR}/${src}_chr${chr}_tagged.vcf.gz"

            # Strip genotypes and compress
            bcftools view -G "$INPUT_VCF" -Oz -o "$SITE_ONLY"
            tabix -f "$SITE_ONLY"

            # Add SOURCE tag inline
            gunzip -c "$SITE_ONLY" | \
                sed "s/IMPUTED;/IMPUTED;SOURCE=${src};/" | \
                bgzip -c > "$MODIFIED_SITE"
            tabix -f "$MODIFIED_SITE"

            # Add header for SOURCE tag
            bcftools annotate \
                -h <(echo '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source imputation panel">') \
                "$MODIFIED_SITE" -Oz -o "$TAGGED_VCF"
            tabix -f "$TAGGED_VCF"

            ANNOTATED_VCFS+=("$TAGGED_VCF")
        else
            echo "⚠️ Warning: missing $INPUT_VCF"
        fi
    done

    # Step 2: Concatenate all site-only VCFs
    CONCAT_VCF="${TMP_DIR}/chr${chr}_concat.vcf.gz"
    bcftools concat -a -Oz -o "$CONCAT_VCF" "${ANNOTATED_VCFS[@]}"
    tabix -f "$CONCAT_VCF"

    # Step 3: Pick best-R² per site
    BEST_SITES="${TMP_DIR}/chr${chr}.best_sites.vcf.gz"
    bcftools view "$CONCAT_VCF" | \
    awk '
    BEGIN {FS="\t"; OFS="\t"}
    /^#/ {print; next}
    {
        key = $1":"$2":"$4":"$5
        match($8, /R2=([^;]+)/, arr)
        r2 = (arr[1] != "") ? arr[1] + 0 : 0
        if (r2 < 0.3) next
        if (!(key in best) || r2 > best_r2[key]) {
            best[key] = $0
            best_r2[key] = r2
        }
    }
    END {
        for (k in best) print best[k]
    }' | bcftools sort -Oz -o "$BEST_SITES"
    tabix -f "$BEST_SITES"

    # Step 4: Extract genotypes from correct original panel using SOURCE filter
    GENOTYPED_VCFS=()
    for src in "${INPUT_DIRS[@]}"; do
        ORIG_VCF="${src}/chr${chr}.dose.vcf.gz"
        SRC_SITES="${TMP_DIR}/${src}_chr${chr}_best_sites.vcf.gz"

        # Filter best_sites for this panel
        bcftools view -i "INFO/SOURCE == \"${src}\"" "$BEST_SITES" -Oz -o "$SRC_SITES"
        tabix -f "$SRC_SITES"

        if [[ -s "$SRC_SITES" && -f "$ORIG_VCF" ]]; then
            OUT_VCF="${TMP_DIR}/${src}_chr${chr}_geno.vcf.gz"
            bcftools view -R "$SRC_SITES" "$ORIG_VCF" -Oz -o "$OUT_VCF"
            tabix -f "$OUT_VCF"
            GENOTYPED_VCFS+=("$OUT_VCF")
        fi
    done

    # Step 5: Merge all genotype files
    MERGED_VCF="${TMP_DIR}/chr${chr}_merged_geno.vcf.gz"
    bcftools merge -Oz -o "$MERGED_VCF" "${GENOTYPED_VCFS[@]}"
    tabix -f "$MERGED_VCF"

    # Step 6: Annotate merged genotypes with SOURCE tag
    FINAL_VCF="${OUTPUT_DIR}/chr${chr}.best_imputed.vcf.gz"
    bcftools annotate -a "$BEST_SITES" -c CHROM,POS,ID,REF,ALT,INFO/SOURCE \
        "$MERGED_VCF" -Oz -o "$FINAL_VCF"
    tabix -f "$FINAL_VCF"

    echo "Done: $FINAL_VCF"
done

# Optional: merge all chromosomes
bcftools concat -a ${OUTPUT_DIR}/chr*.best_imputed.vcf.gz -Oz -o imputed_snps_combined.vcf.gz
tabix -f imputed_snps_combined.vcf.gz

echo "Final output: imputed_snps_combined.vcf.gz"
