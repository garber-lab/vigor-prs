# shuo.shan@umassmed.edu
# this scirpt preps a hg38 vcf.gz file for genotype imputation on Michigan Imputation Server (MIS)
# 1. liftover hg38->hg19; remove variants that failed liftover; 
# 2. correct variants using Will Rayner perl script per MIS advice, for 1000G and HRC reference panels respectively

#!/bin/bash
#BSUB -J pre_imputation_filt
#BSUB -W 8:00
#BSUB -n 1
#BSUB -R "rusage[mem=20000]"
#BSUB -o pre_imputation_filt.out
#BSUB -e pre_imputation_filt.err
#BSUB -q long

# Load necessary modules
module load plink
module load plink2/alpha6.1amd
module load htslib
module load bcftools
module load samtools

# Ensure output directory exists
dir=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/MichiganImputationServer_reimpute
mkdir -p ${dir}/data

# ======= inputs =======
CHAIN="${dir}/data/hg38ToHg19.over.chain.gz"    # UCSC chain
LIFTOVER="${dir}/scripts/liftOver"        # UCSC liftOver binary
HG19_FA="${dir}/data/human_g1k_v37.fasta" # hg19/GRCh37 fasta (no 'chr' in names)
# =================================

# ======= create plink file set of hg38 ======
cd ${dir}/data
vcf=/pi/manuel.garber-umw/human/skin/eQTLs/donor_genotyping/hg38/merged_lowGPreset.vcf.gz
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o tmp.setid.vcf.gz ${vcf} # some SNPs don't have rsIDs, rename SNPs with this method
bcftools index -t tmp.setid.vcf.gz
# Make a 3-column sex file: FID  IID  SEX   (1=male, 2=female, 0=unknown)
bcftools query -l tmp.setid.vcf.gz \
  | awk 'BEGIN{OFS="\t"} {print $1, $1, 1}' > sex.update

plink2 --vcf tmp.setid.vcf.gz --update-sex sex.update --split-par b38 --make-bed --out cohort.hg38.withIDs
BFILE="${dir}/data/cohort.hg38.withIDs"
# =================================

# ======= download hg19 fasta and liftover tools =======
cd ${dir}/data
# Download (GRCh37 / b37 without 'chr' prefix)
wget -O human_g1k_v37.fasta.gz \
  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
REF37=${dir}/data/human_g1k_v37.fasta

# Index for tools
samtools faidx ${REF37}

# Create sequence dictionary (.dict) with GATK inside the container
bsub -q interactive -Is singularity exec /share/pkg/containers/gatk/4.6.0.0/gatk-4.6.0.0.sif \
  gatk CreateSequenceDictionary -R $REF37 -O ${REF37%.fasta}.dict

# get hg38tohg19 liftover chain
wget -O hg38ToHg19.over.chain.gz \
  'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'

# get UCSC liftover binary
cd ${dir}/scripts
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/liftOver
chmod +x liftOver
# =================================

# ======= liftover hg38 --> hg19 =======
# 0) Make a PLINK set with chromosome labels PLINK 1.9 will accept (map PAR1/PAR2/XY -> 23)
cd ${dir}/data
OUTPAR="cohort.hg38.parfix"
cp "${BFILE}.bed" "${OUTPAR}.bed"
cp "${BFILE}.fam" "${OUTPAR}.fam"

awk 'BEGIN{OFS="\t"}
{
  c=$1
  if (c=="PAR1" || c=="PAR2" || c=="XY") c=23;   # treat PAR records as X (coords are on X)
  else if (c=="X")  c=23;
  else if (c=="Y")  c=24;
  else if (c=="MT" || c=="M") c=26;
  $1=c; print
}' "${BFILE}.bim" > "${OUTPAR}.bim"
# 0) Make a BED from the cleaned BIM (BED is 0-based; XY already folded into X above)
awk 'BEGIN{OFS="\t"}
{
  c=$1
  if(c==23 || c=="X")        c="chrX";
  else if(c==24 || c=="Y")   c="chrY";
  else if(c==26 || c=="MT")  c="chrM";
  else                       c="chr" c;
  print c, $4-1, $4, $2;     # chrom, start, end, SNPID
}' "${OUTPAR}.bim" > sites.hg38.bed

# 1) LiftOver to hg19
${LIFTOVER} sites.hg38.bed "$CHAIN" sites.hg19.bed sites.unmapped.bed

# 2) Build PLINK update files from the lifted BED
#    (update-map: ID newBP) (update-chr: ID newCHRcode)
awk 'BEGIN{OFS="\t"}
{
  gsub(/^chr/,"",$1);
  print $4, $3             # ID, new BP
}' sites.hg19.bed > update.bp

#    update-chr: ID newCHRcode
awk 'BEGIN{OFS="\t"} !/^#/ {
  gsub(/^chr/,"",$1)
  chr=$1
  if(chr=="X") c=23; else if(chr=="Y") c=24; else if(chr=="M"||chr=="MT") c=26; else c=chr+0
  print $4, c
}' sites.hg19.bed > update.chr

#    Variants that failed liftover → exclude list
cat sites.unmapped.bed | grep -v "Deleted" | awk '{print $4}' | sort -u > failLift.ids


# 3) Update PLINK fileset to hg19 coords (and drop unmapped)
# Re-split PARs for hg19 (so PAR sites become XY again)
#    (PLINK lets you combine --update-chr and --update-map in one run; must use --make-bed.)
plink --bfile "${OUTPAR}" \
  --exclude failLift.ids \
  --update-map update.bp \
  --update-chr update.chr \
  --update-sex sex.update \
  --make-bed --out cohort.hg19.tmp

#plink2 --bfile cohort.hg19.tmp \
#  --split-par b37 \
#  --snps-only just-acgt \
#  --remove <(awk '$2=="skineQTL-F086" || $2=="skineQTL-F104"{print $1, $2}' cohort.hg19.tmp.fam) \
#  --recode vcf bgz \
#  --output-chr MT \
#  --out impute_ready.hg19

tabix -f -p vcf impute_ready.hg19.vcf.gz
# =================================
# ========= Rayner correct variants (strand/allele/ID/position fixes) ==============
# download Rayner's tools
cd ${dir}/scripts
wget https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
unzip HRC-1000G-check-bim-v4.3.0.zip

cd ${dir}/data
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz # panel sites files
wget https://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz

# remove F086/F104 from hg19 PLINK set (individuals failed plinkQC due to heterozygosity test)
cd ${dir}/data
awk '$2=="skineQTL-F086" || $2=="skineQTL-F104"{print $1,$2}' cohort.hg19.tmp.fam > remove.ids
plink --bfile cohort.hg19.tmp --hwe .001 --remove remove.ids --make-bed --out panelprep
plink --bfile panelprep --freq --out panelprep

# ========================================
# run Rayner against HRC ref panel
mkdir -p ${dir}/pre_imputation_vcf_HRC
cd ${dir}/pre_imputation_vcf_HRC
ln -s ${dir}/data/panelprep.bed .
ln -s ${dir}/data/panelprep.bim .
ln -s ${dir}/data/panelprep.fam .
ln -s ${dir}/data/panelprep.frq .
ln -s ${dir}/data/panelprep.log .
ln -s ${dir}/data/panelprep.nosex .
REF_HRC="${dir}/data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
# running this generates the Run-plink.sh file, which contains the output file name. 
perl ${dir}/scripts/HRC-1000G-check-bim.pl \
    -b panelprep.bim \
    -f panelprep.frq \
    -r "$REF_HRC" \
    -o ${dir}/pre_imputation_vcf_HRC \
    -h
# Apply the fixes Rayner generated
bash Run-plink.sh
# (This creates a corrected PLINK dataset; the final --out name is inside Run-plink.sh.)
mkdir -p vcf
mv *.vcf vcf
cd vcf
for f in *.vcf;do
        echo ${f}
        bcftools view --threads 8 -Oz -o ${f}.gz ${f}
        tabix -f -p vcf ${f}.gz
done
rm *.vcf

# ========================================
# run Rayner against 1000G ref panel
mkdir -p ${dir}/pre_imputation_vcf_1000G
cd ${dir}/pre_imputation_vcf_1000G
ln -s ${dir}/data/panelprep.bed .
ln -s ${dir}/data/panelprep.bim .
ln -s ${dir}/data/panelprep.fam .
ln -s ${dir}/data/panelprep.frq .
ln -s ${dir}/data/panelprep.log .
ln -s ${dir}/data/panelprep.nosex .
REF_1KG="${dir}/data/1000GP_Phase3_combined.legend.gz"
perl ${dir}/scripts/HRC-1000G-check-bim.pl \
    -b ${dir}/data/panelprep.bim \
    -f ${dir}/data/panelprep.frq \
    -r "$REF_1KG" \
    -o ${dir}/pre_imputation_vcf_1000G \
    -g -p ALL
# Apply the fixes Rayner generated
bash Run-plink.sh
# (This creates a corrected PLINK dataset; the final --out name is inside Run-plink.sh.)
mkdir -p vcf
mv *.vcf vcf
cd vcf
for f in *.vcf;do
        echo ${f}
        bcftools view --threads 8 -Oz -o ${f}.gz ${f}
        tabix -f -p vcf ${f}.gz
done
rm *.vcf
# =================================
