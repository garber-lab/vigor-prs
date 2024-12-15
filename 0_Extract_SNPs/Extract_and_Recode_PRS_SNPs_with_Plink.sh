#!/bin/bash

# Run PLINK command
plink \
  --bfile "/pi/manuel.garber-umw/human/VIGOR/tj/ibd/all_samples_qc" \
  --extract "/home/genevieve.roberts-umw/vigor-prs/0_Extract_SNPs/PLINK_extract_SNPs.tsv" \
  --recode A \
  --out "Roberts19_Risk_Score_SNPs_with_Tags"
