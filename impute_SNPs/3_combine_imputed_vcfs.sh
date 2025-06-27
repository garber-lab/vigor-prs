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

    for dir in "${INPUT_DIRS[@]}"; do
        INPUT_VCF="${dir}/chr${chr}.dose.vcf.gz"

        if [[ -f "$INPUT_VCF" ]]; then
            SITE_ONLY="${TMP_DIR}/${dir}_chr${chr}_site.vcf.gz"
            MODIFIED_SITE="${TMP_DIR}/${dir}_chr${chr}_site_modified.vcf.gz"
            TAGGED_VCF="${TMP_DIR}/${dir}_chr${chr}_tagged.vcf.gz"

            # Step 1: Remove genotypes and compress (bgzip)
            bcftools view -G "$INPUT_VCF" -Oz -o "$SITE_ONLY"
            tabix -f "$SITE_ONLY"

            # Step 2: Decompress, inject SOURCE tag, recompress
            gunzip -c "$SITE_ONLY" | sed "s/IMPUTED;/IMPUTED;SOURCE=${dir};/" | bgzip -c > "$MODIFIED_SITE"
            tabix -f "$MODIFIED_SITE"

            # Step 3: Add SOURCE to header and write final annotated VCF
            bcftools annotate \
                -h <(echo '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source imputation panel">') \
                "$MODIFIED_SITE" \
                -Oz -o "$TAGGED_VCF"
            tabix -f "$TAGGED_VCF"

            ANNOTATED_VCFS+=("$TAGGED_VCF")
        else
            echo "Warning: $INPUT_VCF not found"
        fi
    done

    # Step 4: Concatenate all source-tagged VCFs
    CONCAT_VCF="${TMP_DIR}/chr${chr}_concat.vcf.gz"
    bcftools concat -a -Oz -o "$CONCAT_VCF" "${ANNOTATED_VCFS[@]}"
    tabix -f "$CONCAT_VCF"

    # Step 5: Select best R2 per site using awk + bgzip
    bcftools view "$CONCAT_VCF" | \
    awk '
    BEGIN {FS="\t"; OFS="\t"}
    /^#/ {print; next}
    {
        key = $1":"$2":"$4":"$5
        match($8, /R2=([^;]+)/, arr)
        r2 = (arr[1] != "") ? arr[1] + 0 : 0
        if (!(key in best) || r2 > best_r2[key]) {
            best[key] = $0
            best_r2[key] = r2
        }
    }
    END {
        for (k in best) print best[k]
    }' | bcftools sort -Oz -o "${OUTPUT_DIR}/chr${chr}.best_imputed.vcf.gz"


    tabix -f "${OUTPUT_DIR}/chr${chr}.best_imputed.vcf.gz"
done
