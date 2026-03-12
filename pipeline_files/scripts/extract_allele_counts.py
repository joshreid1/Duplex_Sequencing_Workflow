import subprocess
import sys
from collections import Counter
import os
import re
import gzip
import csv
import pysam


def get_sample_id(bam_file):
    """Extract sample ID from BAM file path."""
    basename = os.path.basename(bam_file)
    sample_id = basename.split('.')[0]
    return sample_id


def run_samtools_view(bam_file, tag_value, simplex_file, duplex_file):
    """Run samtools view to filter reads by XW tag and create simplex and duplex BAMs."""
    try:
        subprocess.run(
            ["samtools", "view", "-@4", "-h", "-b", "-d", f"XW:{tag_value}", bam_file,
             "--output-unselected", duplex_file, "-o", simplex_file],
            check=True
        )
        print(f"Created filtered BAMs: {duplex_file} and {simplex_file}")

        subprocess.run(["samtools", "index", duplex_file], check=True)
        subprocess.run(["samtools", "index", simplex_file], check=True)
        print(f"Indexed filtered BAMs: {duplex_file} and {simplex_file}")

    except subprocess.CalledProcessError as e:
        print(f"Error processing BAM file: {e.stderr}", file=sys.stderr)
        sys.exit(1)


def count_fragments_at_position(bam_file, chrom, pos, ref_allele, alt_allele, min_bq=10, min_mapq=1):
    """Count unique fragments supporting ref/alt at a position using pysam."""
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        
        fragment_votes = {}
        
        ref_allele = ref_allele.upper()
        alt_allele = alt_allele.upper()
        
        is_snv = len(ref_allele) == 1 and len(alt_allele) == 1
        is_insertion = len(ref_allele) < len(alt_allele)
        is_deletion = len(ref_allele) > len(alt_allele)
        
        fetch_end = pos + max(len(ref_allele), len(alt_allele))
        for read in bam.fetch(chrom, pos-1, fetch_end):
            if (read.is_unmapped or read.is_duplicate or 
                read.is_secondary or read.is_supplementary):
                continue
            
            if read.mapping_quality < min_mapq:
                continue
            
            aligned_pairs = read.get_aligned_pairs()
            
            allele = None
            base_quality = 0
            
            if is_snv:
                for query_pos, ref_pos in aligned_pairs:
                    if ref_pos == pos - 1 and query_pos is not None:
                        base = read.query_sequence[query_pos].upper()
                        base_qual = read.query_qualities[query_pos]
                        
                        if base_qual >= min_bq:
                            allele = base
                            base_quality = base_qual
                        break
            
            elif is_insertion:
                insertion_seq = alt_allele[len(ref_allele):]
                
                for i, (query_pos, ref_pos) in enumerate(aligned_pairs):
                    if ref_pos == pos - 1 and query_pos is not None:
                        query_bases = []
                        j = i + 1
                        while j < len(aligned_pairs):
                            next_qpos, next_rpos = aligned_pairs[j]
                            if next_rpos is None and next_qpos is not None:
                                query_bases.append(read.query_sequence[next_qpos].upper())
                                j += 1
                            else:
                                break
                        
                        anchor_qual = read.query_qualities[query_pos]
                        ins_quals = [read.query_qualities[query_pos]] + \
                                   [read.query_qualities[aligned_pairs[k][0]] 
                                    for k in range(i+1, i+1+len(query_bases))]
                        
                        if min(ins_quals) >= min_bq:
                            if len(query_bases) == len(insertion_seq) and \
                               ''.join(query_bases) == insertion_seq:
                                allele = alt_allele
                                base_quality = min(ins_quals)
                            else:
                                allele = ref_allele
                                base_quality = anchor_qual
                        break
            
            elif is_deletion:
                deletion_length = len(ref_allele) - len(alt_allele)
                
                for i, (query_pos, ref_pos) in enumerate(aligned_pairs):
                    if ref_pos == pos - 1:
                        del_count = 0
                        j = i + 1
                        while j < len(aligned_pairs):
                            next_qpos, next_rpos = aligned_pairs[j]
                            if next_qpos is None and next_rpos is not None:
                                if next_rpos < pos - 1 + len(ref_allele):
                                    del_count += 1
                                j += 1
                            else:
                                break
                        
                        if query_pos is not None:
                            anchor_qual = read.query_qualities[query_pos]
                            if anchor_qual >= min_bq:
                                if del_count == deletion_length:
                                    allele = alt_allele
                                    base_quality = anchor_qual
                                else:
                                    allele = ref_allele
                                    base_quality = anchor_qual
                        break
            
            if allele is None:
                continue
            
            read_name = read.query_name
            
            if read_name not in fragment_votes:
                fragment_votes[read_name] = (allele, base_quality)
            else:
                existing_allele, existing_qual = fragment_votes[read_name]
                
                if existing_allele == allele:
                    if base_quality > existing_qual:
                        fragment_votes[read_name] = (allele, base_quality)
                else:
                    if base_quality > existing_qual:
                        fragment_votes[read_name] = (allele, base_quality)
        
        bam.close()
        
        allele_counts = Counter()
        for allele, _ in fragment_votes.values():
            allele_counts[allele] += 1
        
        ref_count = allele_counts.get(ref_allele, 0)
        alt_count = allele_counts.get(alt_allele, 0)
        total_coverage = sum(allele_counts.values())
        other_count = total_coverage - ref_count - alt_count
        
        return chrom, pos, ref_count, alt_count, other_count, total_coverage
        
    except Exception as e:
        print(f"Error processing position {chrom}:{pos} in {bam_file}: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return None


def process_vcf(vcf_file):
    """Parse the VCF file and yield variant information."""
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    with open_func(vcf_file, mode) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4].split(',')[0]
            
            ad_info = None
            if len(parts) >= 10:
                format_keys = parts[8].split(':')
                sample_values = parts[9].split(':')
                if "AD" in format_keys:
                    ad_index = format_keys.index("AD")
                    if ad_index < len(sample_values):
                        ad_info = sample_values[ad_index]
            yield chrom, pos, ref, alt, ad_info


def main(bam_file, vcf_file):
    sample_id = get_sample_id(bam_file)
    simplex_bam = f"{sample_id}_simplex.bam"
    duplex_bam = f"{sample_id}_duplex.bam"
    
    run_samtools_view(bam_file, "0", simplex_bam, duplex_bam)

    with open("vaf_info.csv", "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        header = [
            "CHR", "POS", "REF", "ALT", "VCF_AD",
            "BAM_REF", "BAM_ALT", "BAM_OTHER",
            "DUPLEX_REF", "DUPLEX_ALT", "DUPLEX_OTHER",
            "SIMPLEX_REF", "SIMPLEX_ALT", "SIMPLEX_OTHER"
        ]
        csvwriter.writerow(header)

        for chrom, pos, ref, alt, ad_info in process_vcf(vcf_file):
            original_result = count_fragments_at_position(bam_file, chrom, int(pos), ref, alt, 
                                                         min_bq=10, min_mapq=1)
            duplex_result   = count_fragments_at_position(duplex_bam, chrom, int(pos), ref, alt,
                                                         min_bq=10, min_mapq=1)
            simplex_result  = count_fragments_at_position(simplex_bam, chrom, int(pos), ref, alt,
                                                         min_bq=10, min_mapq=1)

            if not original_result or not duplex_result or not simplex_result:
                continue

            _, _, bam_ref, bam_alt, bam_other, _ = original_result
            _, _, duplex_ref, duplex_alt, duplex_other, _ = duplex_result
            _, _, simplex_ref, simplex_alt, simplex_other, _ = simplex_result

            row = [
                chrom, pos, ref, alt, ad_info if ad_info else ".",
                bam_ref, bam_alt, bam_other,
                duplex_ref, duplex_alt, duplex_other,
                simplex_ref, simplex_alt, simplex_other
            ]
            csvwriter.writerow(row)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_allele_counts.py <BAM_FILE> <VCF_FILE>")
        sys.exit(1)

    bam_file, vcf_file = sys.argv[1:3]
    main(bam_file, vcf_file)