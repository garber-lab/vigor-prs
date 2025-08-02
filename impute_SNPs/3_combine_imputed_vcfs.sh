#!/bin/bash
#BSUB -J combine_imputed_vcfs[1-22]
#BSUB -R "rusage[mem=64000]"
#BSUB -o combine_imputed_vcfs_%I.out
#BSUB -e combine_imputed_vcfs_%I.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

# Load necessary modules
module load plink2/alpha6.1amd
module load bcftools
module load htslib

cd /home/genevieve.roberts-umw/imputed_genotypes

mkdir -p logs
OUTPUT_DIR="./combined_output"
mkdir -p "$OUTPUT_DIR"

TMP_DIR="./imputation_tmp"
mkdir -p "$TMP_DIR"

INPUT_DIRS=("1KG" "HRC" "MHC_alleles")
MASTER_FAM="/pi/manuel.garber-umw/human/VIGOR/tj/ibd/run_2/all_samples_ibd.fam"

CHROM_LIST=({1..22})
chr="${CHROM_LIST[$((LSB_JOBINDEX-1))]}"

echo "Processing chromosome $chr..."

PGEN_BASENAMES=()
R2_SUMMARY_FILES=()

# Step 1: Convert VCFs to plink2 pgen format with conditional allele renaming
for src in "${INPUT_DIRS[@]}"; do
    INPUT_VCF="${src}/chr${chr}.dose.vcf.gz"
    if [[ -f "$INPUT_VCF" ]]; then
        OUT_PREFIX="${src}_chr${chr}"

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

        if [[ "$src" == "MHC_alleles" && "$chr" == "6" ]]; then
            echo "Handling MHC_alleles chr6 with HLA_ variant exception..."

            HLA_VCF="${TMP_DIR}/${src}_chr${chr}_HLA_only.vcf.gz"
            bcftools view -i 'ID ~ "^HLA_"' "$INPUT_VCF" -O z -o "$HLA_VCF"
            bcftools index -f "$HLA_VCF"

            plink2 --vcf "$HLA_VCF" \
                --extract-if-info "R2 >= 0.3" \
                --new-id-max-allele-len 10 missing \
                --update-ids "${TMP_DIR}/${src}_chr${chr}_samples.txt" \
                --make-pgen --out "${TMP_DIR}/${src}_chr${chr}"
        else
            plink2 --vcf "$INPUT_VCF" \
            --extract-if-info "R2 >= 0.3" \
            --set-all-var-ids '@:#:$r:$a' \
            --new-id-max-allele-len 10 missing \
            --update-ids "${TMP_DIR}/${src}_chr${chr}_samples.txt" \
            --make-pgen --out "${TMP_DIR}/${OUT_PREFIX}"
        fi
            
        PGEN_BASENAMES+=("${TMP_DIR}/${OUT_PREFIX}")

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
mapfile -t R2_SUMMARY_FILES < <(find "${TMP_DIR}" -name '*_variant_r2.tsv')

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

# Step 3: Intersect samples and create filtered BFILES
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
        
        plink2 --pfile "${TMP_DIR}/${OUT_PREFIX}" \
            --extract "$VARIANT_LIST" \
            --make-just-pvar --out "${TMP_DIR}/${OUT_PREFIX}_filtered"

        FILTERED_BFILES+=("${TMP_DIR}/${OUT_PREFIX}_filtered")
    else
        echo "Warning: missing ${TMP_DIR}/${OUT_PREFIX}.pgen"
    fi
done

# Step 4: Merge BFILES
MERGE_LIST="${TMP_DIR}/chr${chr}_merge_list.txt"
> "$MERGE_LIST"
for f in "${FILTERED_BFILES[@]:1}"; do
    echo "$f" >> "$MERGE_LIST"
done

MERGED_PREFIX="${TMP_DIR}/chr${chr}_merged"

plink --bfile "${FILTERED_BFILES[0]}" \
    --merge-list "$MERGE_LIST" \
    --make-bed \
    --out "$MERGED_PREFIX"

plink2 --bfile "$MERGED_PREFIX" \
    --make-just-pvar \
    --out "$MERGED_PREFIX"

# Step 5: Annotate .pvar file
FILTERED_PVAR_FILES=()
for f in "${FILTERED_BFILES[@]}"; do
    FILTERED_PVAR_FILES+=("${f}.pvar")
done

awk -v r2_file="${BEST_R2_FILE}" '
    BEGIN {
        FS = OFS = "\t"
        while ((getline < r2_file) > 0) {
            id = $1
            r2 = $2
            src = $3
            info[id] = "R2=" r2 ";SOURCE=" src
        }
    }
    !/^#/ {
        id = $3
        if (id in info) {
            print $3, info[id]
        }
    }
' "${FILTERED_PVAR_FILES[@]}" | sort -k1,1V -k2,2n | uniq > "${MERGED_PREFIX}.filtered_combined.anno"

MERGED_FAM="${MERGED_PREFIX}.fam"
UPDATED_SEX_FILE="${TMP_DIR}/chr${chr}_sex_update.txt"
UPDATED_PHENO_FILE="${TMP_DIR}/chr${chr}_pheno_update.txt"

awk '{print $1, $2}' "$MERGED_FAM" > "${TMP_DIR}/chr${chr}_samples_to_update.txt"

grep -Ff "${TMP_DIR}/chr${chr}_samples_to_update.txt" "$MASTER_FAM" | awk '{print $1, $2, $5}' > "$UPDATED_SEX_FILE"
grep -Ff "${TMP_DIR}/chr${chr}_samples_to_update.txt" "$MASTER_FAM" | awk '{print $1, $2, $6}' > "$UPDATED_PHENO_FILE"

plink2 \
--bfile "$MERGED_PREFIX" \
--update-sex "$UPDATED_SEX_FILE" \
--pheno "$UPDATED_PHENO_FILE" \
--make-pgen \
--out "${OUTPUT_DIR}/chr${chr}_imputed_combined"

awk '
BEGIN { FS = OFS = "\t" }
NR==FNR { anno[$1] = $2; next }
/^#/ {
    if ($1 == "#CHROM") {
        print $1, $2, $3, $4, $5, "INFO"
    } else {
        print
    }
    next
}
{
    info = ($3 in anno) ? anno[$3] : "R2=NA;SOURCE=NA"
    print $1, $2, $3, $4, $5, info
}
' "${MERGED_PREFIX}.filtered_combined.anno" "${MERGED_PREFIX}.pvar" > "${OUTPUT_DIR}/chr${chr}_imputed_combined.pvar"

# Add header lines
pvar_file="${OUTPUT_DIR}/chr${chr}_imputed_combined.pvar"
{
    echo '##INFO=<ID=R2,Number=1,Type=Float,Description="Imputation INFO score from reference panel">'
    echo '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Imputation reference panel: HRC, 1KG, MHC_alleles">'
    cat "$pvar_file"
} > tmp && mv tmp "$pvar_file"

# Step 6: Create and post-process VCF
FINAL_PLINK_OUT_PREFIX="${TMP_DIR}/chr${chr}_imputed_combined"
plink2 --pfile "${OUTPUT_DIR}/chr${chr}_imputed_combined" \
    --export vcf bgz \
    --out "$FINAL_PLINK_OUT_PREFIX"

CLEANED_VCF="${TMP_DIR}/tmp_chr${chr}_imputed_combined_cleaned.vcf.gz"
FINAL_VCF="${OUTPUT_DIR}/chr${chr}_imputed_combined_hg19.vcf.gz"

bcftools view -e 'REF ~ "^<"' "${FINAL_PLINK_OUT_PREFIX}.vcf.gz" -O z -o "$CLEANED_VCF"
bcftools index -f "$CLEANED_VCF"

bcftools annotate \
--rename-chrs <(echo -e "${chr}\tchr${chr}") \
-O z -o "$FINAL_VCF" "$CLEANED_VCF"
bcftools index -f "$FINAL_VCF"
