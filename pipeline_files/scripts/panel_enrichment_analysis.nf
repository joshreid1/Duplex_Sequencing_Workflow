#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline parameters
 */
params.controls = null      // TSV file with control sample IDs (1 per row, no header)
params.cases = null         // TSV file with case sample IDs (1 per row, no header)
params.gene = null          // Gene symbol (e.g., KMT2D)
params.vcf_dir = '/vast/scratch/users/reid.j/austin_panel/'
params.outdir = 'results'

/*
 * Validate required parameters
 */
if (!params.controls || !params.cases || !params.gene) {
    error "Missing required parameters. Please provide: --controls, --cases, --gene"
}

/*
 * Gene coordinate lookup (hg38)
 * Expand this map with additional genes as needed
 */
def getGeneCoordinates(gene) {
    def coords = [
        'KMT2D': 'chr12:49018978-49060794',
        'BRCA1': 'chr17:43044295-43125483',
        'TP53': 'chr17:7668421-7687490',
        'SCN1A': 'chr2:166845706-166930414',
        'MECP2': 'chrX:154021573-154137105',
		'MTOR': 'chr1:11106535-11262551',
        // Add more genes here
    ]
    
    if (!coords.containsKey(gene)) {
        error "Gene ${gene} not found in coordinate lookup. Please add coordinates manually."
    }
    
    return coords[gene]
}

/*
 * Process: Locate VCF files for samples
 */
process LOCATE_VCF_FILES {
    tag "${group}"
    
    input:
    tuple val(group), path(sample_list)
    
    output:
    tuple val(group), path("${group}_vcf_list.txt")
    
    script:
    """
    #!/bin/bash
    
    # Read sample IDs and find corresponding VCF files
    while IFS= read -r sample_id; do
        # Skip empty lines
        [[ -z "\$sample_id" ]] && continue
     
        # Search for VCF file recursively, selecting the most recent
        vcf_file=\$(find ${params.vcf_dir} -type f -name "\${sample_id}*.hard-filtered.vcf.gz" -exec ls -t {} + | head -n 1)

        if [[ -n "\$vcf_file" ]]; then
            # Check if index file exists
            tbi_file="\${vcf_file}.tbi"
            if [[ ! -f "\$tbi_file" ]]; then
                echo "WARNING: Index file not found for \${vcf_file}" >&2
            fi
            echo "\${sample_id}\t\${vcf_file}\t\${tbi_file}" >> ${group}_vcf_list.txt
        else
            echo "WARNING: VCF file not found for sample \${sample_id}" >&2
        fi
    done < ${sample_list}
    
    # Check if any VCF files were found
    if [[ ! -f ${group}_vcf_list.txt ]]; then
        echo "ERROR: No VCF files found for ${group} samples" >&2
        exit 1
    fi
    """
}


/*
 * Process: Subset VCF by gene region and count variants
 */
process SUBSET_AND_COUNT_VARIANTS {
    tag "${sample_id}"
    
    input:
    tuple val(group), val(sample_id), path(vcf), path(tbi), val(region)
    
    output:
    tuple val(group), val(sample_id), path("${sample_id}.${params.gene}.vcf.gz"), path("${sample_id}.${params.gene}.vcf.gz.tbi"), path("${sample_id}.variant_count.txt")
    
    script:
    """
    # Subset VCF to gene region
    bcftools view -r ${region} ${vcf} | bcftools view -e 'FORMAT/AF > 0.30' | bcftools view -i 'FORMAT/AD[0:1] > 2' -Oz -o ${sample_id}.${params.gene}.vcf.gz  
    bcftools index -t ${sample_id}.${params.gene}.vcf.gz
    
    # Count variants (exclude reference genotypes)
    variant_count=\$(bcftools view -H ${sample_id}.${params.gene}.vcf.gz | wc -l)
    echo -e "${sample_id}\t\${variant_count}" > ${sample_id}.variant_count.txt
    """
}


/*
 * Process: Enrichment analysis (Fisher's exact test)
 */
process ENRICHMENT_ANALYSIS {
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple path(control_counts), path(case_counts), val(gene), val(region)
    
    output:
    path "*enrichment_results.txt"
    path "*enrichment_summary.txt"
    
    script:
    """
    #!/usr/bin/env python3
    import scipy.stats as stats
    import numpy as np
    
    # Read counts
    def read_counts(filename):
        samples_with_variants = 0
        samples_without_variants = 0
        total_variants = 0
        
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split('\\t')
                if len(parts) != 2:
                    continue
                sample, count = parts
                count = int(count)
                total_variants += count
                if count > 0:
                    samples_with_variants += 1
                else:
                    samples_without_variants += 1
        
        return samples_with_variants, samples_without_variants, total_variants
    
    # Read control and case counts
    ctrl_with, ctrl_without, ctrl_total_vars = read_counts('${control_counts}')
    case_with, case_without, case_total_vars = read_counts('${case_counts}')
    
    # Total samples
    ctrl_total = ctrl_with + ctrl_without
    case_total = case_with + case_without
    
    # Create contingency table
    # Rows: [cases, controls]
    # Columns: [with variants, without variants]
    table = [[case_with, case_without],
             [ctrl_with, ctrl_without]]
    
    # Perform Fisher's exact test
    oddsratio, pvalue = stats.fisher_exact(table, alternative='greater')
    
    # Calculate proportions
    case_proportion = case_with / case_total if case_total > 0 else 0
    ctrl_proportion = ctrl_with / ctrl_total if ctrl_total > 0 else 0
    
    # Write detailed results
    with open('${gene}_enrichment_results.txt', 'w') as f:
        f.write("=== Enrichment Analysis Results ===\\n")
        f.write(f"Gene: ${gene}\\n")
        f.write(f"Region: ${region}\\n\\n")
        
        f.write("--- Cases ---\\n")
        f.write(f"Total samples: {case_total}\\n")
        f.write(f"Samples with variants: {case_with}\\n")
        f.write(f"Samples without variants: {case_without}\\n")
        f.write(f"Total variants: {case_total_vars}\\n")
        f.write(f"Proportion with variants: {case_proportion:.4f}\\n\\n")
        
        f.write("--- Controls ---\\n")
        f.write(f"Total samples: {ctrl_total}\\n")
        f.write(f"Samples with variants: {ctrl_with}\\n")
        f.write(f"Samples without variants: {ctrl_without}\\n")
        f.write(f"Total variants: {ctrl_total_vars}\\n")
        f.write(f"Proportion with variants: {ctrl_proportion:.4f}\\n\\n")
        
        f.write("--- Statistical Test ---\\n")
        f.write(f"Fisher's Exact Test (one-sided, greater)\\n")
        f.write(f"Odds Ratio: {oddsratio:.4f}\\n")
        f.write(f"P-value: {pvalue:.6e}\\n\\n")
        
        f.write("--- Interpretation ---\\n")
        if pvalue < 0.05:
            f.write(f"Variants in ${gene} are significantly enriched in cases (p < 0.05)\\n")
        else:
            f.write(f"No significant enrichment detected (p >= 0.05)\\n")
    
    # Write summary
    with open('${gene}_enrichment_summary.txt', 'w') as f:
        f.write(f"Gene\\tRegion\\tCases_Total\\tCases_With_Var\\tControls_Total\\tControls_With_Var\\tOdds_Ratio\\tP_value\\tSignificant\\n")
        sig = "Yes" if pvalue < 0.05 else "No"
        f.write(f"${gene}\\t${region}\\t{case_total}\\t{case_with}\\t{ctrl_total}\\t{ctrl_with}\\t{oddsratio:.4f}\\t{pvalue:.6e}\\t{sig}\\n")
    """
}

/*
 * Workflow
 */
workflow {
    // Get gene coordinates
    gene_region = getGeneCoordinates(params.gene)
    
    // Create channels from input files and add group labels
    controls_ch = Channel.fromPath(params.controls)
        .map { tuple('controls', it) }
    
    cases_ch = Channel.fromPath(params.cases)
        .map { tuple('cases', it) }
    
    // Combine both groups
    all_samples = controls_ch.mix(cases_ch)
    
    // Locate VCF files
    vcf_lists = LOCATE_VCF_FILES(all_samples)
    
    // Parse VCF lists and create sample channels with VCF and index files
    sample_vcfs = vcf_lists
        .map { group, vcf_list ->
            vcf_list.text.readLines().collect { line ->
                def parts = line.split('\t')
                def sample_id = parts[0]
                def vcf_path = parts[1]
                def tbi_path = parts[2]
                tuple(group, sample_id, file(vcf_path), file(tbi_path), gene_region)
            }
        }
        .flatMap()
        .view { group, sample_id, vcf_file, tbi_file, gene_region ->
        "Sample: ${sample_id} | VCF: ${vcf_file}"
    }
    
    // Subset and count variants for each sample
    variant_counts = SUBSET_AND_COUNT_VARIANTS(sample_vcfs)
    
    // Separate control and case counts and aggregate
    control_counts_list = variant_counts
        .filter { it[0] == 'controls' }
        .map { it[4] }  // Get the count file (5th element in tuple)
        .collectFile(name: 'controls_counts.txt', newLine: false, sort: true)
    
    case_counts_list = variant_counts
        .filter { it[0] == 'cases' }
        .map { it[4] }  // Get the count file (5th element in tuple)
        .collectFile(name: 'cases_counts.txt', newLine: false, sort: true)
    
    // Combine for enrichment analysis
    enrichment_input = control_counts_list
        .combine(case_counts_list)
        .map { ctrl, cases -> tuple(ctrl, cases, params.gene, gene_region) }
    
    // Run enrichment analysis
    ENRICHMENT_ANALYSIS(enrichment_input)
}

workflow.onComplete {
    println "Pipeline completed!"
    println "Results available in: ${params.outdir}"
}
