#!/usr/bin/env python3
"""
Script to calculate depth of coverage at every nucleotide position in gene regions
from a BAM file, restricted to reads with MAPQ > 0.
"""

import pysam
import pandas as pd
import argparse
import os
from pathlib import Path

def parse_bed_file(bed_file):
    """
    Parse BED file and return list of regions.
    
    Args:
        bed_file (str): Path to BED file
        
    Returns:
        list: List of tuples (chrom, start, end)
    """
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chrom = fields[0]
                start = int(fields[1])  # BED is 0-based
                end = int(fields[2])    # BED end is exclusive
                regions.append((chrom, start, end))
    return regions

def calculate_coverage(bam_file, regions, min_mapq=1):
    """
    Calculate coverage depth at each position in the specified regions.
    
    Args:
        bam_file (str): Path to BAM file
        regions (list): List of tuples (chrom, start, end)
        min_mapq (int): Minimum mapping quality threshold
        
    Returns:
        list: List of tuples (chrom, pos, depth)
    """
    coverage_data = []
    
    # Open BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        
        for chrom, start, end in regions:
            print(f"Processing region {chrom}:{start}-{end}")
            
            # Get pileup for the region
            for pileupcolumn in bamfile.pileup(chrom, start, end, 
                                             truncate=True, 
                                             stepper='samtools',
                                             ignore_overlaps=True,
                                             ignore_orphans=True,
                                             min_base_quality= min_mapq):
                
                # Filter reads by mapping quality
                depth = 0
                for pileupread in pileupcolumn.pileups:
                    if (not pileupread.is_del and 
                        not pileupread.is_refskip):
                        depth += 1
                
                # Convert to 1-based coordinates for output
                pos = pileupcolumn.pos + 1
                coverage_data.append((chrom, pos, depth))
    
    return coverage_data

def extract_sample_id(bam_file):
    """
    Extract sample ID from BAM filename.
    
    Args:
        bam_file (str): Path to BAM file
        
    Returns:
        str: Sample ID extracted from filename
    """
    filename = Path(bam_file).stem
    # Remove .bam extension if present
    if filename.endswith('.bam'):
        filename = filename[:-4]
    return filename

def main():
    parser = argparse.ArgumentParser(description='Calculate coverage depth from BAM file for gene regions')
    parser.add_argument('--bam', required=True, help='Path to BAM file')
    parser.add_argument('--bed', required=True, help='Path to BED file with gene regions')
    parser.add_argument('--ref', help='Path to reference FASTA file (optional)')
    parser.add_argument('--min-mapq', type=int, default=1, help='Minimum mapping quality (default: 1)')
    parser.add_argument('--output', help='Output file path (default: auto-generated)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.bam):
        raise FileNotFoundError(f"BAM file not found: {args.bam}")
    if not os.path.exists(args.bed):
        raise FileNotFoundError(f"BED file not found: {args.bed}")
    if args.ref and not os.path.exists(args.ref):
        raise FileNotFoundError(f"Reference FASTA file not found: {args.ref}")
    
    # Check if BAM index exists
    bam_index = args.bam + '.bai'
    if not os.path.exists(bam_index):
        print(f"Creating BAM index: {bam_index}")
        pysam.index(args.bam)
    
    # Extract sample ID from filename
    sample_id = extract_sample_id(args.bam)
    print(f"Sample ID: {sample_id}")
    
    # Parse BED file
    print("Parsing BED file...")
    regions = parse_bed_file(args.bed)
    print(f"Found {len(regions)} regions")
    
    # Calculate coverage
    print("Calculating coverage...")
    coverage_data = calculate_coverage(args.bam, regions, args.min_mapq)
    
    # Create DataFrame
    df = pd.DataFrame(coverage_data, columns=['chrom', 'pos', 'depth'])
    
    # Generate output filename if not provided
    if args.output:
        output_file = args.output
    else:
        output_file = f"{sample_id}_coverage_depth.tsv"
    
    # Save to file
    print(f"Saving results to: {output_file}")
    df.to_csv(output_file, sep='\t', index=False, header=False)
    
    # Print summary statistics
    total_positions = len(df)
    covered_positions = len(df[df['depth'] > 0])
    mean_depth = df['depth'].mean()
    median_depth = df['depth'].median()
    
    print(f"\nSummary Statistics:")
    print(f"Total positions: {total_positions:,}")
    print(f"Covered positions: {covered_positions:,} ({covered_positions/total_positions*100:.2f}%)")
    print(f"Mean depth: {mean_depth:.2f}")
    print(f"Median depth: {median_depth:.2f}")
    print(f"Max depth: {df['depth'].max()}")

if __name__ == "__main__":
    main()
