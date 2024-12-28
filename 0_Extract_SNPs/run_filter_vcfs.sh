#!/bin/bash
#BSUB -J filter_vcfs[1-140]
#BSUB -R rusage[mem=4096]
#BSUB -o filter_vcfs_%J_%I.out
#BSUB -e filter_vcfs_%J_%I.err
#BSUB -q short
#BSUB -n 1

# Load necessary modules
module load bcftools

# Define input directory and temporary directory
INPUT_DIR="/pi/manuel.garber-umw/human/VIGOR"
TMP_DIR="/pi/manuel.garber-umw/human/VIGOR/groberts/filtered_vcfs"

# Create temporary directory for filtered VCFs
mkdir -p $TMP_DIR
chmod -R u+rwx $TMP_DIR

# Get list of VCF files
VCF_FILES=($(find $INPUT_DIR -name "*.vcf.gz"))

# Get the current VCF file for this job
vcf_file=${VCF_FILES[$((LSB_JOBINDEX-1))]}
sample_id=$(basename $vcf_file | sed 's/\.vcf\.gz//')

# Filter VCF by INFO/RAF > 0.001
output_vcf="$TMP_DIR/${sample_id}_filtered.vcf.gz"
echo "Processing $vcf_file for $sample_id"
bcftools view -i 'INFO/RAF > 0.001' $vcf_file -Oz -o $output_vcf
chmod u+rwx $output_vcf
bcftools index $output_vcf
