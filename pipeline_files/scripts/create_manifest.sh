#!/bin/bash

# Define the output manifest file
OUTPUT_FILE=manifest.tsv

# Check if directories were provided as arguments
if [ "$#" -eq 0 ]; then
    echo "Please provide one or more directories to search."
    exit 1
fi

# Write the header to the manifest file
echo -e "ID\tBATCH\tVCF\tBAM" > "$OUTPUT_FILE"

# Loop through each directory provided as argument
for SEARCH_DIR in "$@"; do
    if [ ! -d "$SEARCH_DIR" ]; then
        echo "Directory $SEARCH_DIR does not exist, skipping."
        continue
    fi

    # Use find command to locate all .hard-filtered.vcf.gz files and process them
    find "$SEARCH_DIR" -name "*.hard-filtered.vcf.gz" | while read -r filepath; do
        # Extract ID and BATCH information from the file path
        ID=$(basename "$filepath" .hard-filtered.vcf.gz)
        BATCH=$(echo "$filepath" | grep -oE '(MBC|EST)[0-9]+' | awk '!a[$0]++' | paste -sd'-')

        echo $filepath
        echo $ID
        echo $BATCH
        
        # Search in the same directory for a .bam file that contains the ID in its filename
        BAM=$(find "$(dirname "$filepath")" -maxdepth 1 -type f -name "*${ID}*.bam" | head -n 1)
        
        # Append the entry to the manifest file with BAM as the fourth column
        echo -e "${ID}\t${BATCH}\t${filepath}\t${BAM}" >> "$OUTPUT_FILE"
    done
done

echo "Manifest file generated: $OUTPUT_FILE"
