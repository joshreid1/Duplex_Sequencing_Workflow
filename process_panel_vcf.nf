#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Input parameters
params.sample_info			= "${projectDir}/pipeline_files/manifest_files/test_sample_info.tsv"
params.bed_file             = "${projectDir}/pipeline_files/gene_lists/austin_panel_targets.bed"
params.ref_fasta			= '/stornext/Bioinf/data/lab_bahlo/projects/epilepsy/hg38/reference/fasta/Homo_sapiens_assembly38.fasta'
params.spliceai_distance    = 500

// Duplex/simplex read support thresholds for filtering somatic variants.
//   A variant passes if it meets ANY of the three conditions:
//   1. Strong duplex support alone:         DUPLEX_ALT >= min_duplex
//   2. Combined duplex + simplex support:   DUPLEX_ALT >= min_duplex_with_simplex & SIMPLEX_ALT >= min_simplex_with_duplex
//   3. Strong simplex support alone:        SIMPLEX_ALT >= min_simplex

params.min_duplex              = 2
params.min_simplex             = 3
params.min_duplex_with_simplex = 1
params.min_simplex_with_duplex = 2

// vcfanno toml files
params.clinvar_toml         		= "${projectDir}/pipeline_files/vcfanno_files/clinvar_20250330.toml"
params.cosmic_toml          		= "${projectDir}/pipeline_files/vcfanno_files/cosmic_20251203.toml"
params.gnomad_toml        			= "${projectDir}/pipeline_files/vcfanno_files/gnomad_v4.0.0.toml"
params.gnomad_postprocess_toml 		= "${projectDir}/pipeline_files/vcfanno_files/gnomad_v4.0.0_postprocess.toml"
params.spliceai_lua					= "${projectDir}/pipeline_files/vcfanno_files/spliceai.lua"
params.spliceai_postprocess_toml 	= "${projectDir}/pipeline_files/vcfanno_files/spliceai_postprocess.toml"


// vep files
params.vep_cache_dir		= "/vast/projects/bahlo_cache/vep_cache/"
params.vep_alphamissense 	= '/vast/projects/bahlo_cache/annotation/alphamissense/AlphaMissense_hg38.tsv.gz'
params.vep_revel			= '/vast/projects/bahlo_cache/annotation/REVEL/revel_1.3.hg38.vep.tsv.gz'

process gnomad {

	tag "${ID} ${GROUP}"

    container 'quay.io/biocontainers/vcfanno:0.2.6--0'

	cpus = 1
	memory = { 1 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	input:
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	output:
		path ('gnomad.vcf')
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
			
	script:
	"""
	#gnomAD annotations 
	vcfanno -p 4 ${params.gnomad_toml} ${VCF} > tmp.vcf
	vcfanno -p 4 ${params.gnomad_postprocess_toml} tmp.vcf > gnomad.vcf
	rm tmp.vcf
	"""
} 	

process Pre_Filter_Variants {

	tag "${ID} ${GROUP}"

    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

	cpus = 1
	memory = { 1 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	input:
		path (vcf)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
	output:
		path ('filtered_variants.vcf')
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
			
	script:
	"""
	bcftools view -i 'FILTER=="PASS"' ${vcf} | \
	bcftools view -e 'AF_gnomad_total_4.0 > 0.001' -Ov -o filtered_variants.vcf 
	#bcftools view -e 'FORMAT/AF > 0.25'
	"""
}

process VEP {

	tag "${ID} ${GROUP}"

    container = 'quay.io/biocontainers/ensembl-vep:115.2--pl5321h2a3209d_1'

	cpus = 1
	memory = { 1 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	input:
		path (vcf)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
	output:
		path ('vep.vcf')
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

    script:
	"""
	vep --cache --dir ${params.vep_cache_dir} --cache_version 115 --assembly GRCh38 \
			-i ${vcf} -o vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs \
			--fasta ${params.ref_fasta} --offline --sift b --polyphen b --ccds --hgvs --hgvsg --symbol \
			--numbers --protein --variant_class --pick_allele_gene --force_overwrite \
			--plugin AlphaMissense,file=${params.vep_alphamissense},transcript_match=1 \
			--plugin REVEL,file=${params.vep_revel}
    """
}

process CADD_Run_Container {
	tag "${ID} ${GROUP}"

	container 'oras://docker.io/joshreid1/cadd-scoring:v1.6_edit'
	containerOptions '-B /vast/projects/bahlo_epilepsy/somatic_annotation_data/CADD-scripts/data/annotations:/CADD-scripts/data/annotations --writable-tmpfs'

	cpus = 1
	memory = { 10 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}
		
	input:
		path (vcf)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
	
		
	output:
		path(vcf)
		path("tmp.tsv.gz")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
	
	shell:
	'''
	#Link local software
	ln -s /stornext/System/data/software/rhel/9/base/tools/snakemake/8.11.3/bin/snakemake /usr/local/bin/snakemake
	ln -s /stornext/System/data/software/rhel/9/base/bioinf/bcftools/1.20/bin/bcftools /usr/local/bin/bcftools

	if [[ $(bcftools query -f '%ALT\n' !{vcf} | uniq) == "*" ]]; then
		cp !{vcf} tmp.vcf
		continue
	else
		grep -v "^#" !{vcf} | cut -f 1-5 | sed 's/^chr//' > tmp.vcf
		/CADD-scripts/CADD.sh tmp.vcf
	fi
	'''
}


process Process_CADD {
	tag "${ID} ${GROUP}"

    container 'community.wave.seqera.io/library/htslib_vcfanno:8044b99f5458cd69'

	cpus = 1
	memory = { 10 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	input:
		path(vcf)
		path(cadd_tsv) 
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
		
	output:
		tuple path("pass.annotated.cadd.vcf.gz"), path("pass.annotated.cadd.vcf.gz.tbi")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
	
	shell:
	'''
	cp "$(readlink -f !{cadd_tsv})" ./cadd_tmp

	tabix -f -b 2 -e 2 -s 1 ./cadd_tmp

	cadd_path=$(realpath ./cadd_tmp)
	echo '[[annotation]]' > cadd.toml
	echo "file= '${cadd_path}'" >> cadd.toml
	echo 'names=["CADD_Score"]' >> cadd.toml
	echo 'ops=["mean"]' >> cadd.toml
	echo 'columns=[6]' >> cadd.toml

	vcfanno cadd.toml !{vcf} > pass.annotated.cadd.vcf

	bgzip pass.annotated.cadd.vcf
	tabix pass.annotated.cadd.vcf.gz  
	'''
}
 
process ClinVar {
	tag "${ID} ${GROUP}"

    container 'quay.io/biocontainers/vcfanno:0.2.6--0'	

	input:
		tuple path(vcf), path(vcf_index)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
				
	output:
		path("*.clinvar.vcf")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	shell:
	'''
	vcfanno  !{params.clinvar_toml} !{vcf} > !{ID}_!{GROUP}.clinvar.vcf
	'''
}

process Cosmic {
	tag "${ID} ${GROUP}"

    container 'quay.io/biocontainers/vcfanno:0.2.6--0'		

	input:
		path(vcf)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
				
	output:
		path("*.cosmic.vcf")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	shell:
	'''
	vcfanno  !{params.cosmic_toml} !{vcf} > !{ID}_!{GROUP}.cosmic.vcf
	'''
}

process SpliceAI_Run {
	tag "${ID} ${GROUP}"	

	container 'community.wave.seqera.io/library/python_pip_keras_setuptools_pruned:1c71801b2a7b49db'

	cpus = 1
	memory = { 16 * task.attempt + ' GB' }
	time = { 2 * task.attempt + ' h'}

	input:
		path(vcf)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	output:
		path("*.spliceai.vcf")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	script:
	""" 
	export PYTHON_EGG_CACHE=./
	spliceai -I ${vcf} \
	-O ${ID}_${GROUP}.spliceai.vcf \
	-R ${params.ref_fasta} \
	-D ${params.spliceai_distance} \
	-A /vast/projects/bahlo_epilepsy/somatic_annotation_data/SpliceAI/gencode.v43.canonical.annotation.txt
	"""
}

process Process_SpliceAI {
	tag "${ID} ${GROUP}"

    container 'quay.io/biocontainers/vcfanno:0.2.6--0'

	cpus = 1
	memory = { 2 * task.attempt + ' GB' }
	time = { 2 * task.attempt + ' h'}

	input:
		path(vcf)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	output:
		path("*_processed_spliceai.vcf")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	script:
	""" 
	vcfanno -lua ${params.spliceai_lua} ${params.spliceai_postprocess_toml} ${vcf} \
		> ${ID}_${GROUP}_processed_spliceai.vcf
	"""
}

process Filter_Variants_Vembrane {
	
    tag "${ID} ${GROUP}"

    container 'quay.io/biocontainers/vembrane:2.5.0--pyhdfd78af_0'
    
    publishDir "filtered_variants", mode: 'copy'
    
    input:
    path(vcf)
    tuple val(GROUP), val(ID), val(VCF), val(BAM)
    
    output:
    path("*_candidate_variants.vcf")
    tuple val(GROUP), val(ID), val(VCF), val(BAM)

    script:
	"""

	# Vembrane filter logic:
	# INCLUDES variants that are NOT benign/likely benign in ClinVar AND meet at least ONE of:
	#   a) High-impact consequence (stop gain/loss, start loss, frameshift, inframe indels)
	#   b) Pathogenic/Likely pathogenic in ClinVar
	#   c) CADD score > 20
	#   d) SpliceAI score > 0.5
	#   e) REVEL score > 0.6

	cat > vembrane_expr.txt <<'EXPR'
(not "benign" in str(INFO.get("ClinVarSIG", "")).lower() and not "benign" in str(INFO.get("ClinVarSIGCONF", "")).lower() and (any(any(cons in ann.split("|")[1] for cons in ["stop_gained", "stop_lost", "start_lost", "frameshift_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant"]) for ann in INFO.get("CSQ", [])) or "pathogenic" in str(INFO.get("ClinVarSIG", "")).lower() or ((INFO.get("CADD_Score") or 0) > 20) or (any(max([float(s) if s and s != "." else 0 for s in spliceai.split("|")[2:6]]) > 0.5 for spliceai in INFO.get("SpliceAI", []) if len(spliceai.split("|")) > 5)) or any(float(str(ANN["REVEL"])) > 0.6 for _ in [1] if str(ANN["REVEL"]) not in ("", ".", "NoValue()"))))
EXPR

    vembrane filter \
		--annotation-key CSQ \
        --overwrite-number-info "COSMIC_Sample_Count=." \
        --output ${ID}_candidate_variants.vcf \
        "\$(cat vembrane_expr.txt)" \
        ${vcf}
    """
}

process Compress_Index {

	tag "${ID} ${GROUP}"	

	cpus = 1
	memory = { 1 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	input:
		path(candidate_variants)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	output:
		tuple path("*.vcf.gz"), path("*.vcf.gz.tbi")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	script:
	"""
	bgzip ${candidate_variants}
	tabix ${candidate_variants}.gz
	"""
}

process Check_Simplex_Duplex_VAF {

	tag "${ID} ${GROUP}"
	
	cpus = 1
	memory = { 16 * task.attempt + ' GB' }
	time = { 4 * task.attempt + ' h'}

	container 'community.wave.seqera.io/library/pysam_samtools:e22dac7b8f8c32bf'

	// Publish only the BAM files and their indexes
	//publishDir "updated_results", mode: "copy", pattern: "*.{bam,bai}"

	input:
		tuple path(candidate_variants_gz), path(candidate_variants_tbi)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	output:
		tuple path("*simplex.bam"), path("*simplex.bam.bai"), emit: simplex_bam
		tuple path("*duplex.bam"),  path("*duplex.bam.bai"),  emit: duplex_bam
		tuple path("*candidate_variants.vcf.gz", includeInputs:true), path("*candidate_variants.vcf.gz.tbi", includeInputs:true)
		path("vaf_info.csv")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	script:
	"""
	python ${projectDir}/pipeline_files/scripts/extract_allele_counts.py ${BAM} ${candidate_variants_gz}
	"""
}

process Filter_Variants_VAF {
	cpus = 1
	memory = { 16 * task.attempt + ' GB' }
	time = { 4 * task.attempt + ' h'}

	container 'community.wave.seqera.io/library/bioconductor-variantannotation_r-tidyverse:bc8ea2e5386b79c1'

	input:
		tuple path(simplex_bam), path(simplex_index)
		tuple path(duplex_bam),  path(duplex_index)
		tuple path(candidate_variants_gz), path(candidate_variants_tbi)
		path(vaf_info_csv)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	output:
		path("*.RDS"), emit: rds
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	script:
	"""
	Rscript ${projectDir}/pipeline_files/scripts/filter_austin_panel_variants.R \
		${candidate_variants_gz} ${vaf_info_csv} ${ID}_${GROUP}.RDS ${projectDir} \
		${params.min_duplex} ${params.min_simplex} \
		${params.min_duplex_with_simplex} ${params.min_simplex_with_duplex}
	"""
}

process Generate_Report {
	publishDir "filtered_variants", mode: "copy"

	memory = '32 GB'

	container 'community.wave.seqera.io/library/r-dplyr_r-openxlsx_zip:04bb05dca48c0236'

	input:
	path rds_files

	output:
	path "*.xlsx"

	script:
	"""
	Rscript ${projectDir}/pipeline_files/scripts/combine_austin_panel_rds.R ./ filtered_variants.xlsx
	"""
}

process Check_Coverage {
	cpus = 2
	memory = '4 GB'

	input: tuple val(GROUP), val(ID), val(VCF), val(BAM)

	output: tuple val(GROUP), val(ID), val(VCF), val(BAM), emit: meta
			path ("*.gene-list.per-base.bed"), emit: bed

	tag "${ID} ${GROUP} coverage"

	script:
	"""
	python ${projectDir}/pipeline_files/scripts/check_coverage.py --bam ${BAM} --bed ${params.bed_file} --output ${ID}.${GROUP}.gene-list.per-base.bed
	"""
}

process Plot_Depths {
	cpus = 2
	memory = { 10 * task.attempt + ' GB' }

	publishDir "coverage_data", mode: 'copy'

	input:
	file bed_files

	output:
	path "*.png"
	path "coverage_summary.tsv"

	script:
	"""
	Rscript ${projectDir}/pipeline_files/scripts/plot_panel_depths.R .
	"""
}


workflow {
    samples_ch = Channel
        .fromPath(params.sample_info)
        .splitCsv(header: true, sep: "\t")
        .map { row ->
            [
                Group        : row.Group,
                Sample_ID    : row.Sample_ID,
                VCF_Filepath : row.VCF_Filepath,
                BAM_Filepath : row.BAM_Filepath
            ]
        }

    samples_ch
        .collect()
        .view { rows ->
            def total = rows.size()
            def groupCounts = rows.countBy { it.Group }

            def summary = groupCounts
                .collect { group, n -> "${n} ${group} samples" }
                .join(', ')

            """Processing ${total} samples
${summary}"""
        }

    sampleChannel = samples_ch.map { row ->
        tuple(row.Group, row.Sample_ID, row.VCF_Filepath, row.BAM_Filepath)
    }

	// Variant annotation and filtering
	gnomad(sampleChannel)
	| Pre_Filter_Variants
	| VEP
	| CADD_Run_Container
	| Process_CADD
	| ClinVar
	| Cosmic
	| SpliceAI_Run
	| Process_SpliceAI
	| Filter_Variants_Vembrane
	| Compress_Index
	| Check_Simplex_Duplex_VAF
	| Filter_Variants_VAF
	Filter_Variants_VAF.out.rds.collect().set { all_rds_files }
	Generate_Report(all_rds_files)
	
	// Coverage analysis and plotting
	Check_Coverage(sampleChannel)
	Plot_Depths(Check_Coverage.out.bed.collect())

}