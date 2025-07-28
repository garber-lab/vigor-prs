#!/bin/bash
#BSUB -J liftover_chr22
#BSUB -R "rusage[mem=48000]"
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
OUT_VCF_WITH_CHR=${OUT_DIR}/${VCF_PREFIX}_addchr.dose.vcf.gz
REJECT_VCF=${OUT_DIR}/test_hg38_chr${CHR}_rejected_variants.vcf
LIFTED_VCF=${OUT_DIR}/test_hg38_chr${CHR}.dose.vcf

# ---- STEP 1: Add chr prefixes to contigs ----
echo "Adding 'chr' prefixes to contigs for chr${CHR}..."
bcftools annotate \
  --rename-chrs <(echo -e "${CHR}\tchr${CHR}") \
  -O z -o "$OUT_VCF_WITH_CHR" "$IN_VCF"

bcftools index "$OUT_VCF_WITH_CHR"

# ---- STEP 2: Run Picard LiftoverVcf ----
echo "Running LiftoverVcf for chr${CHR}..."
JAVA_MEM=40g
java -Xmx$JAVA_MEM -jar $PICARDJAR LiftoverVcf \
    I="$OUT_VCF_WITH_CHR" \
    O="$LIFTED_VCF" \
    CHAIN="$CHAIN" \
    REJECT="$REJECT_VCF" \
    R="$REFERENCE" \
    MAX_RECORDS_IN_RAM=1000000

echo "Liftover complete for chr${CHR}."

#module load samtools/1.3

## Get the chain file
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
#gunzip hg19ToHg38.over.chain.gz

## Get the reference genome
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#gunzip hg38.fa.gz

## Make the dictionary and index for the hg38 reference genome
#java -jar $PICARDJAR CreateSequenceDictionary R=hg38.fa O=hg38.dict
# samtools faidx hg38.fa
