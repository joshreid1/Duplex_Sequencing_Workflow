#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Input parameters
params.sample_info			= '/vast/scratch/users/reid.j/duplex_sequencing_workflow/mbc002_ids_plus_vcfpath_plus_bam.tsv'
params.bed_file             = "${projectDir}/pipeline_files/gene_lists/austin_panel_targets.bed"
params.spliceai_distance    = 500

// vcfanno toml files
params.clinvar_toml         		= "${projectDir}/pipeline_files/vcfanno_files/clinvar_20250330.toml"
params.cosmic_toml          		= "${projectDir}/pipeline_files/vcfanno_files/cosmic_20251203.toml"
params.gnomad_toml        			= "${projectDir}/pipeline_files/vcfanno_files/gnomad_v4.0.0.toml"
params.gnomad_postprocess_toml 		= "${projectDir}/pipeline_files/vcfanno_files/gnomad_v4.0.0_postprocess.toml"
params.spliceai_lua					= "${projectDir}/pipeline_files/vcfanno_files/spliceai.lua"
params.spliceai_postprocess_toml 	= "${projectDir}/pipeline_files/vcfanno_files/spliceai_postprocess.toml"

params.ref_fasta			= '/stornext/Bioinf/data/lab_bahlo/projects/epilepsy/hg38/reference/fasta/Homo_sapiens_assembly38.fasta'
params.mosdepth				= '/stornext/Bioinf/data/lab_bahlo/software/apps/mosdepth/mosdepth_v0.3.8' 
params.vep_cache_dir 		= "/stornext/Bioinf/data/lab_bahlo/ref_db/vep-cache/"

// Import modules
//include { SpliceAI_Run } from '/stornext/Bioinf/data/lab_bahlo/users/reid.j/Splice_Pipeline/modules.nf'

process gnomad {

	tag "${ID} ${GROUP}"

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

	module 'bcftools/1.20'

	tag "${ID} ${GROUP}"

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
	module 'ensembl-vep/112'

	tag "${ID} ${GROUP}"

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
	vep --cache --dir ${params.vep_cache_dir} --cache_version 104 --assembly GRCh38 \
			-i ${vcf} -o vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs \
			--fasta ${params.ref_fasta} --offline --sift b --polyphen b --ccds --hgvs --hgvsg --symbol \
			--numbers --protein --af --af_1kg --max_af --variant_class --pick_allele_gene --force_overwrite
	"""
}

process CADD_Run_Container {
	tag "${ID} ${GROUP}"

	cpus = 1
	memory = { 10 * task.attempt + ' GB' }
	time = { 1 * task.attempt + ' h'}

	container '/vast/projects/reidj-project/containers/cadd-scoring_latest.sif'
	//container 'community.wave.seqera.io/library/cadd-scripts:1.7.3--e283e32d163015db'

	containerOptions '-B /stornext/Bioinf/data/lab_bahlo/users/reid.j/cadd/CADD-scripts-master/data/annotations:/CADD-scripts/data/annotations --writable-tmpfs'
		
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

	input:
		tuple path(vcf), path(vcf_index)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
				
	output:
		tuple path("*.clinvar.vcf.gz"), path("*.clinvar.vcf.gz.tbi")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	shell:
	'''
	vcfanno  !{params.clinvar_toml} !{vcf} > !{ID}_!{GROUP}.clinvar.vcf
	bgzip !{ID}_!{GROUP}.clinvar.vcf
	tabix !{ID}_!{GROUP}.clinvar.vcf.gz
	'''
}

process Cosmic {
	tag "${ID} ${GROUP}"	

	input:
		tuple path(vcf), path(vcf_index)
		tuple val(GROUP), val(ID), val(VCF), val(BAM)
				
	output:
		tuple path("*.cosmic.vcf.gz"), path("*.cosmic.vcf.gz.tbi")
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	shell:
	'''
	vcfanno  !{params.cosmic_toml} !{vcf} > !{ID}_!{GROUP}.cosmic.vcf

	bgzip !{ID}_!{GROUP}.cosmic.vcf
	tabix !{ID}_!{GROUP}.cosmic.vcf.gz
	'''
}

process SpliceAI_Run {
	tag "${ID} ${GROUP}"	

	cpus = 1
	memory = { 16 * task.attempt + ' GB' }
	time = { 2 * task.attempt + ' h'}

	container '/stornext/Bioinf/data/lab_bahlo/users/reid.j/analysis/Apptainer/spliceai.img'

	input:
		tuple path(vcf), path(vcf_index)
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
	-A /stornext/Bioinf/data/lab_bahlo/users/reid.j/spliceai/gencode.v43.canonical.annotation.txt
	"""
}

process Process_SpliceAI {
	tag "${ID} ${GROUP}"

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
	vcfanno -lua ${params.spliceai_lua} ${params.spliceai_postprocess_toml} ${vcf} > ${ID}_${GROUP}_processed_spliceai.vcf
	"""
}

process Filter_Variants_Vembrane {
	
    tag "${ID} ${GROUP}"
    
    label 'vembrane'

    container 'quay.io/biocontainers/vembrane:1.0.7--pyhdfd78af_0'
    
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
	#   b) Pathogenic/likely pathogenic in ClinVar
	#   c) CADD score >20
	#   d) SpliceAI score >0.8
	# TO ADD: Revel > 0.6

	cat > vembrane_expr.txt <<'EXPR'
(not "benign" in str(INFO.get("ClinVarSIG", "")).lower() and not "benign" in str(INFO.get("ClinVarSIGCONF", "")).lower() and (any(any(cons in ann.split("|")[1] for cons in ["stop_gained", "stop_lost", "start_lost", "frameshift_variant", "inframe_insertion", "inframe_deletion", "protein_altering_variant"]) for ann in INFO.get("CSQ", [])) or "pathogenic" in str(INFO.get("ClinVarSIG", "")).lower() or ((INFO.get("CADD_Score") or 0) > 20) or (any(max([float(s) if s and s != "." else 0 for s in spliceai.split("|")[2:6]]) > 0.8 for spliceai in INFO.get("SpliceAI", []) if len(spliceai.split("|")) > 5))))
EXPR

    vembrane filter \
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
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	output:
		path("*.RDS"), emit: rds
		tuple val(GROUP), val(ID), val(VCF), val(BAM)

	script:
	"""
	Rscript ${projectDir}/pipeline_files/scripts/filter_austin_panel_variants.R ${candidate_variants_gz} vaf_info.csv ${ID}_${GROUP}.RDS ${projectDir}   
	"""
}

process Generate_Report {
	publishDir "updated_results", mode: "copy"

	memory = '32 GB'

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

	publishDir 'mosdepth_results', mode: 'copy', pattern: "*.bed"

	input: tuple val(GROUP), val(ID), val(VCF), val(BAM)

	output: tuple val(GROUP), val(ID), val(VCF), val(BAM)
			path ("*.gene-list.per-base.bed")

	tag "${ID} ${GROUP} coverage"

	script:
	"""
	python ${projectDir}/pipeline_files/scripts/check_coverage.py --bam ${BAM} --bed ${params.bed_file} --output ${ID}.${GROUP}.gene-list.per-base.bed
	"""
}

workflow {
	// Create a channel from the sample info file
	Channel
		.fromPath(params.sample_info)
		.splitCsv(header: true, sep: "\t")
		.map { row -> [row.Group, row.Sample_ID, row.VCF_Filepath, row.BAM_Filepath] }
		.view()
		.set { sampleChannel }


	// Process the sample channel
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

	Check_Coverage(sampleChannel)

	// Collect all RDS files
	Filter_Variants_VAF.out.rds.collect().set { all_rds_files }

	// Pass the collected RDS files to Generate_Report
	Generate_Report(all_rds_files)
}


/////

process Join_VCF {
	cpus = 2
	memory = { 10 * task.attempt + ' GB' }
	time = { 2 * task.attempt + ' h'}

	module 'bcftools/1.20:htslib'

	input:
		path('vcf', arity: '1..*')

	output:
		tuple path("*.gz"), path("*.tbi")

	shell:
	'''
	if [ $(ls vcf* 2>/dev/null | wc -l) -gt 1 ]; then
		snpsift split -j vcf* | bcftools sort --temp-dir ./ -Oz -o final.vcf.gz --write-index=tbi
	else
		bcftools sort --temp-dir ./ -Oz -o final.vcf.gz --write-index=tbi vcf*
	fi
	'''	
}

process Split_Vcf {
	label 'C1M1T1'
	
	input:
		path(vcf) 

	output:
		path("merged.*.vcf")

	shell:
	'''
	variant_number=`expr $(zgrep -v "^#" !{vcf} | wc -l)`
	
	if [ $variant_number -gt 1000 ]
	then
		task_number=500
	elif [ $variant_number -gt 10 ]
	then
		task_number=10
	else
		task_number=2
	fi
	
	split_number=`expr $variant_number / $task_number`

	split_number=`expr $variant_number / $task_number`
	if [ $split_number -lt 5 ]; then
		split_number=5
	fi

	echo "Variant number = $variant_number" >> error_check.txt
	echo "Number of tasks = $task_number" >> error_check.txt
	echo "Number of variants per split = $split_number" >> error_check.txt

	snpsift split -l $split_number !{vcf} 
	'''
}