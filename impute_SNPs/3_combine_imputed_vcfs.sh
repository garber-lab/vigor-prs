#!/bin/bash
#BSUB -J combine_imp_pgen
#BSUB -R "rusage[mem=8192]"
#BSUB -o combine_imp_pgen.out
#BSUB -e combine_imp_pgen.err
#BSUB -q short
#BSUB -W 1:00
#BSUB -n 1

# Load necessary modules
module load plink2/alpha6.1amd
module load bcftools
module load htslib

cd /home/genevieve.roberts-umw/imputed_genotypes

OUTPUT_DIR="./combined_output"
mkdir -p "$OUTPUT_DIR"

TMP_DIR="./imputation_tmp"
mkdir -p "$TMP_DIR"

INPUT_DIRS=("1KG" "HRC" "MHC_alleles")

for chr in {21..22}; do
    echo "Processing chromosome $chr..."

    PGEN_BASENAMES=()
    R2_SUMMARY_FILES=()

    # Step 1: Convert VCFs to plink2 pgen format with allele renaming
    for src in "${INPUT_DIRS[@]}"; do
        INPUT_VCF="${src}/chr${chr}.dose.vcf.gz"
        if [[ -f "$INPUT_VCF" ]]; then
            OUT_PREFIX="${src}_chr${chr}"

            # Step 1a: Extract sample names and build FID/IID map
            bcftools query -l "$INPUT_VCF" | awk -F'_' '{
                old_fid = 0;
                old_iid = $0;
                new_fid = $1"_"$2;
                new_iid = "";
                for (i = 3; i <= NF; i++) {
                    new_iid = (new_iid == "") ? $i : new_iid"_"$i;
                }
                print old_fid, old_iid, new_fid, new_iid;
            }' > "${TMP_DIR}/${src}_chr${chr}_samples.txt"

            # Step 1b: Convert to PGEN with FID/IID assignment
            plink2 --vcf "$INPUT_VCF" \
            --extract-if-info "R2 >= 0.3" \
            --set-all-var-ids '@:#:$r:$a' \
            --new-id-max-allele-len 10 missing \
            --update-ids "${TMP_DIR}/${src}_chr${chr}_samples.txt" \
            --make-pgen --out "${TMP_DIR}/${OUT_PREFIX}"

            PGEN_BASENAMES+=("${TMP_DIR}/${OUT_PREFIX}")

            # Step 1c: Save variant R2 values for later merging
            awk -v src="$src" '
            BEGIN { OFS = "\t" }
            !/^#/ {
                info = $6
                split(info, fields, ";")
                r2 = -1
                for (i in fields) {
                    if (fields[i] ~ /^R2=/) {
                        split(fields[i], kv, "=")
                        r2 = kv[2] + 0
                    }
                }
                if (r2 >= 0.3 && $3 != ".") {
                    print $3, r2, src
                }
            }' "${TMP_DIR}/${OUT_PREFIX}.pvar" > "${TMP_DIR}/${OUT_PREFIX}_variant_r2.tsv"
        else
            echo "Warning: missing $INPUT_VCF"
        fi
    done

    # Step 2: Find best-R² variant per site
    BEST_R2_FILE="${TMP_DIR}/chr${chr}_best_r2.tsv"

    R2_SUMMARY_FILES=($(ls "${TMP_DIR}"/*_variant_r2.tsv))
    
    cat "${R2_SUMMARY_FILES[@]}" | \
    awk '
    BEGIN { OFS = "\t" }
    {
        key = $1
        if (key in best_r2) {
            if ($2 > best_r2[key]) {
                best_r2[key] = $2
                best_variant[key] = $0
            }
        } else {
            best_r2[key] = $2
            best_variant[key] = $0
        }
    }
    END {
        for (k in best_variant) {
            print best_variant[k]
        }
    }' > "${BEST_R2_FILE}"
    sort "${BEST_R2_FILE}" -o "${BEST_R2_FILE}"

    # Step 3: Extract variant list (plink2 ID format), intersect samples, and make filtered BFILES
    FILTERED_BFILES=()
    INTERSECT_FILE="${TMP_DIR}/chr${chr}_intersect_samples.txt"
    FIRST=true

    for src in "${INPUT_DIRS[@]}"; do
        OUT_PREFIX="${src}_chr${chr}"
        if [[ -f "${TMP_DIR}/${OUT_PREFIX}.pgen" ]]; then
            VARIANT_LIST="${TMP_DIR}/${OUT_PREFIX}_filtered_variants.txt"
            grep "${src}" "${BEST_R2_FILE}" > "$VARIANT_LIST" 

            SAMPLE_FILE="${TMP_DIR}/${OUT_PREFIX}_samples.txt"
            tail -n +2 "${TMP_DIR}/${OUT_PREFIX}.psam" | awk '{print $1, $2}' > "$SAMPLE_FILE"

            if [ "$FIRST" = true ]; then
                cp "$SAMPLE_FILE" "$INTERSECT_FILE"
                FIRST=false
            else
                grep -F -f "$SAMPLE_FILE" "$INTERSECT_FILE" > "${INTERSECT_FILE}.tmp"
                mv "${INTERSECT_FILE}.tmp" "$INTERSECT_FILE"
            fi
        fi
    done

    for src in "${INPUT_DIRS[@]}"; do
        OUT_PREFIX="${src}_chr${chr}"
        if [[ -f "${TMP_DIR}/${OUT_PREFIX}.pgen" ]]; then
            VARIANT_LIST="${TMP_DIR}/${OUT_PREFIX}_filtered_variants.txt"

            plink2 --pfile "${TMP_DIR}/${OUT_PREFIX}" \
                --extract "$VARIANT_LIST" \
                --keep "$INTERSECT_FILE" \
                --make-bed --out "${TMP_DIR}/${OUT_PREFIX}_filtered"

            FILTERED_BFILES+=("${TMP_DIR}/${OUT_PREFIX}_filtered")
        else
            echo "Warning: missing ${TMP_DIR}/${OUT_PREFIX}.pgen"
            continue
        fi
    done

    # Step 4: Merge all filtered BFILES (pgen merging is not fully developed yet, so do this with plink-1.9 .bed/.bim/.fam files)
    MERGE_LIST="${TMP_DIR}/chr${chr}_merge_list.txt"
    > "$MERGE_LIST"
    for f in "${FILTERED_BFILES[@]:1}"; do
        echo "$f" >> "$MERGE_LIST"
    done

    echo "Files to be merged for chr${chr}:"
    cat "$MERGE_LIST"

    MERGED_PREFIX="${OUTPUT_DIR}/chr${chr}_merged"

    plink --bfile "${FILTERED_BFILES[0]}" \
        --merge-list "$MERGE_LIST" \
        --make-bed --out "$MERGED_PREFIX"
done
