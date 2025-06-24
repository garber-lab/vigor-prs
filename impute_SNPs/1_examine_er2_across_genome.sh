#!/bin/bash
#BSUB -J zgrep_ER2
#BSUB -R "rusage[mem=4096]"
#BSUB -o zgrep_ER2.out
#BSUB -e zgrep_ER2.err
#BSUB -q short
#BSUB -n 1

# Set the working directory
cd /home/genevieve.roberts-umw/imputed_genotypes/1KG

OUTPUT_FILE="filtered_ER2_output.gz"

# Initialize an empty temporary file
TEMP_FILE=$(mktemp)

# Add header to the output
echo -e "CHR\tPOS\tSNP\tA1\tA2\tSTART\tSTOP\tTYPED\tIMPUTED\tAF\tMAF\tAVG_CS\tR2\tER2" > "$TEMP_FILE"

# Loop over chromosome files (chr1.info.gz to chr22.info.gz)
for CHR in {1..22}; do
  FILE="chr${CHR}.info.gz"
  
  # Check if the file exists
  if [[ -f "$FILE" ]]; then
    # Use zgrep to find lines containing "ER2", remove lines starting with "##INFO", and format the output
    zgrep "ER2" "$FILE" | awk -F'\t' '!/^##INFO/ && $1 ~ /^[0-9]+$/ {
      split($8, arr, ";");
      typed = arr[1];
      imputed = arr[2];
      af = substr(arr[3], 4);
      maf = substr(arr[4], 5);
      avg_cs = substr(arr[5], 8);
      r2 = substr(arr[6], 4);  # Remove "R2="
      er2 = substr(arr[7], 5); # Remove "ER2="
      print $1, $2, $3, $4, $5, ".", ".", typed, imputed, af, maf, avg_cs, r2, er2;
    }' OFS="\t" >> "$TEMP_FILE"
  else
    echo "File $FILE not found, skipping..."
  fi
done

# Compress the concatenated output
gzip -c "$TEMP_FILE" > "$OUTPUT_FILE"

# Clean up temporary file
rm "$TEMP_FILE"

echo "Filtered output saved to $OUTPUT_FILE"
