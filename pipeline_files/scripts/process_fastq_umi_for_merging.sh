#!/bin/bash
#SBATCH --job-name=fastq
#SBATCH --mem=32GB
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joshua.reid@unimelb.edu.au

# Process all R1 files in the current directory
for r1_file in *_R1.fastq.gz; do
    # Check if file exists (handles case where no matches found)
    [ -f "$r1_file" ] || continue
    
    # Derive R2 filename
    r2_file="${r1_file/_R1.fastq.gz/_R2.fastq.gz}"
    
    # Check if corresponding R2 file exists
    if [ ! -f "$r2_file" ]; then
        echo "Warning: R2 file not found for $r1_file, skipping..."
        continue
    fi
    
    # Get base name (remove _R1_001.fastq.gz)
    base="${r1_file/_R1.fastq.gz/}"
    
    # Define temporary and output filenames
    temp_r1="${base}_R1_temp.fastq.gz"
    temp_r2="${base}_R2_temp.fastq.gz"
    output_r1="${base}_R1.UMI.fastq.gz"
    output_r2="${base}_R2.UMI.fastq.gz"
    
    echo "Processing: $base"
    
    # Run fastp
    fastp \
        -i "$r1_file" \
        -I "$r2_file" \
        -o "$temp_r1" \
        -O "$temp_r2" \
        -U --umi_loc=per_read \
        --umi_len=3 \
        --umi_skip=2 \
        --disable_adapter_trimming \
        --disable_quality_filtering \
        --disable_length_filtering \
        --disable_trim_poly_g \
        -A -Q -L
    
    # Replace underscore with plus sign in UMI field
    echo "Fixing UMI delimiter for $base..."
    zcat "$temp_r1" | \
        sed 's/:\([ACGTN][ACGTN][ACGTN]\)_\([ACGTN][ACGTN][ACGTN]\) /:\1+\2 /' | \
        gzip > "$output_r1"
    
    zcat "$temp_r2" | \
        sed 's/:\([ACGTN][ACGTN][ACGTN]\)_\([ACGTN][ACGTN][ACGTN]\) /:\1+\2 /' | \
        gzip > "$output_r2"
    
    # Remove temporary files
    rm "$temp_r1" "$temp_r2"
    
    echo "Completed: $output_r1 and $output_r2"
    echo ""
done

echo "All files processed!"
