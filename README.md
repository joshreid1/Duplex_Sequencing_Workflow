# Austin Pathology High‑Depth Duplex Panel Sequencing

## Overview

This repository describes the sequencing workflow and downstream processing 
for the **Austin Pathology high‑depth duplex panel**. 
Data is received from Austin Pathology and stored in the Mediaflux project 
`/projects/proj-5430_epilepsy_ont-1128.4.765/data/austin_panel_data`.

## Library preparation & sequencing

> "DNA libraries were prepared using the SureSelect XT HS2 DNA Kit
> (Agilent) with duplex unique molecular identifiers (UMIs) and a custom
> gene panel targeting 140 genes, including 36 Ras/Raf/MAPK pathway
> genes, 12 mTOR pathway genes, 7 other focal epilepsy genes and 85
> additional cancer genes (see Supplementary Table 2 for gene list).
> Sequencing was performed on an Illumina NextSeq 1000 platform at Austin
> Pathology with 150 bp paired‑end reads.  Read alignment to the
> GRCh38 human reference genome and variant calling were performed using
> the DRAGEN DNA Pipeline on Illumina Basespace (v4.3.13).  UMI tags
> were used for read collapsing during alignment, achieving a median
> consensus depth of 582× per sample across panel positions."

Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12668077/#:~:text=High%2DDepth%20Panel%20Sequencing>

### DRAGEN pipeline documentation

<https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/GPipelineIntro_fDG.htm>

## Data Flow

1. Receive hard‑drive from Austin Pathology (Wendi Lin/Gabi Bradshaw) containing:
   - FASTQ files (currently for storage only, as DRAGEN alignment is performed at Austin Pathology)
   - DRAGEN‑processed files (VCFs, BAMs, QC reports, etc.)

2. Transfer data to Mediaflux:
   ```sh
   /projects/proj-5430_epilepsy_ont-1128.4.765/data/austin_panel_data
   ```
   Several access protocols are available – see  
   <https://github.com/joshreid1/Mediaflux_Tips>

# Variant Annotation & Filtering Pipeline

A Nextflow pipeline for annotating, filtering, and prioritising germline & somatic variants. 

---

## Overview

This pipeline takes per-sample VCF and BAM files as input and performs variant annotation followed by tiered filtering to produce a prioritised candidate variant report. Depth of coverage across target regions is also assessed.

```
Input TSV (sample_info)
       │
       ▼
  gnomAD annotation
       │
  Pre-filter (PASS + gnomAD AF < 0.001)
       │
  VEP annotation
       │
  CADD scoring
       │
  ClinVar annotation
       │
  COSMIC annotation
       │
  SpliceAI scoring
       │
  Vembrane candidate variant filtering
       │
  VAF check (simplex and duplex) 
       │
  Excel output
```

Depth of coverage across panel regions (custom pysam script) also run on all samples.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 22.0
- [Singularity/Apptainer](https://apptainer.org/) (for containerised processes)
- HPC cluster with SLURM (or compatible scheduler)

---

## Input

The pipeline requires a **tab-separated sample sheet** with the following columns:

| Column | Description |
|---|---|
| `Group` | Sample group/cohort identifier (e.g. "MTLE", "MND" or "Control") |
| `Sample_ID` | Unique sample identifier |
| `VCF_Filepath` | Absolute path to the input VCF |
| `BAM_Filepath` | Absolute path to the aligned BAM file |

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `params.sample_info` | (required) | Path to input TSV sample sheet |
| `params.bed_file` | `pipeline_files/gene_lists/austin_panel_targets.bed` | BED file of target regions |
| `params.spliceai_distance` | `500` | SpliceAI distance metric |
| `params.ref_fasta` | GRCh38 FASTA | Reference genome |
| `params.vep_cache_dir` | VEP cache directory | Offline VEP cache path |
| `params.clinvar_toml` | ClinVar vcfanno config | ClinVar annotation TOML |
| `params.cosmic_toml` | COSMIC vcfanno config | COSMIC annotation TOML |
| `params.gnomad_toml` | gnomAD v4.0.0 vcfanno config | gnomAD allele frequency TOML |

---

## Pipeline Steps

### 1. `gnomad`
Annotates each VCF with gnomAD v4.0.0 population allele frequencies using [`vcfanno`](https://github.com/brentp/vcfanno), with a post-processing step to compute derived fields.

### 2. `Pre_Filter_Variants`
Filters variants using `bcftools` to retain only:
- Variants with `FILTER == PASS`
- Variants with gnomAD total AF < 0.1% (`AF_gnomad_total_4.0 > 0.001` excluded)

### 3. `VEP`
Runs [Ensembl VEP](https://www.ensembl.org/vep) (v115.2, GRCh38) in offline mode to annotate functional consequences (SO terms), HGVS notation, scores (SIFT, PolyPhen, AlphaMissense, REVEL), allele frequencies, variant class, and canonical transcript selection via --pick_allele_gene.

### 4. `CADD_Run_Container`
Scores variant deleteriousness using [CADD](https://cadd.gs.washington.edu/) inside a Singularity container. Variants where all ALT alleles are `*` (spanning deletions) are skipped.

### 5. `Process_CADD`
Integrates CADD scores back into the VCF using a dynamically generated `vcfanno` TOML config, then compresses and indexes the output with `bgzip`/`tabix`.

### 6. `ClinVar`
Annotates variants with [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) clinical significance fields via `vcfanno`.

### 7. `Cosmic`
Annotates variants with [COSMIC](https://cancer.sanger.ac.uk/cosmic) somatic variant data via `vcfanno`.

### 8. `SpliceAI_Run`
Predicts splice-site-disrupting variants using [SpliceAI](https://github.com/Illumina/SpliceAI) (±500 bp window by default) inside a Singularity container with GENCODE v43 canonical annotations.

### 9. `Process_SpliceAI`
Post-processes SpliceAI output with `vcfanno` and a Lua script to reformat scores into the INFO field.

### 10. `Filter_Variants_Vembrane`
Applies compound Python-expression filtering using [`vembrane`](https://github.com/vembrane/vembrane). A variant is retained as a **candidate** if it is **not** benign/likely benign in ClinVar AND meets at least one of:
- High-impact consequence: `stop_gained`, `stop_lost`, `start_lost`, `frameshift_variant`, `inframe_insertion`, `inframe_deletion`, or `protein_altering_variant`
- ClinVar significance contains "pathogenic"
- CADD score > 20
- SpliceAI max delta score > 0.5
- REVEL score > 0.6

### 11. `Compress_Index`
Compresses filtered VCF with `bgzip` and indexes with `tabix`.

### 12. `Check_VAF`
Extracts fragment-level allele counts at each variant position from a duplex sequencing BAM, stratified by **all reads**, **duplex reads**, and **simplex reads**. 

1. **Split BAM by read type** — Uses `samtools view -d XW:0` to separate simplex (`XW:0`) from duplex (`XW>0`) reads into two indexed sub-BAMs.
2. **Parse VCF** — Iterates variants (supports `.vcf.gz`), extracting `CHROM`, `POS`, `REF`, `ALT` (first allele only), and `FORMAT/AD`.
3. **Count fragments** — For each variant, counts unique fragment support (by read name) in all three BAMs. Handles SNVs, insertions, and deletions. Reads below `MAPQ 1` or base quality 10 are excluded; duplicate/secondary/supplementary reads are skipped. Where paired reads disagree, the higher base-quality read wins.

## 13. `Filter_Variants_VAF`
Parses an annotated VCF and `vaf_info.csv`, integrates all annotation layers, and applies duplex/simplex coverage filters to produce a per-sample `.RDS` file for `Generate_Report`.
### Default:  `DUPLEX_ALT ≥ 2`, OR `DUPLEX_ALT ≥ 1 & SIMPLEX_ALT ≥ 2`, OR `SIMPLEX_ALT ≥ 3`.

### 14. `Generate_Report`
Collects all per-sample `.RDS` files and generates a consolidated Excel report (`.xlsx`) via an R script.

### 15. `Check_Coverage` *(parallel)*
Runs independently on the input BAM for each sample using (via a Python script) to compute per-base coverage over the target BED regions. Results are published to `coverage_data`.

---

## Outputs

| Directory | Contents |
|---|---|
| `filtered_variants/` | Per-sample candidate variant VCFs |
| `updated_results/` | Per-sample BAM/BAI files + final `filtered_variants.xlsx` report |
| `coverage_data/` | Per-sample per-base coverage BED files |

---

## Usage

```bash
nextflow run main.nf \
  --sample_info /path/to/samples.tsv \
  -profile singularity \
  -resume
```

### Test data:

```bash
nextflow run main.nf \
     --sample_info pipeline_files/manifest_files/test_sample_info.tsv \
     -profile singularity \
     -resume
```

---

## Notes

- All processes use dynamic memory/time retry scaling (`task.attempt`).
- Containerised processes (CADD, SpliceAI, Vembrane) use Singularity images.


## Notes

- Alignment would ideally be run in‑house, however DRAGEN is currently best-practice 
  for handling duplex UMIs and is licensed for use at Austin Pathology.
- `nf-core/sarek` only supports single UMIs as of writing.  See issue:
  <https://github.com/nf-core/sarek/issues/1630>.  The Agilent AGeNT tool
  is hoped to appear in a future release.

- Early runs suffered from low coverage.  This was resolved by
  increasing DNA input (100 ng) and limiting to eight samples per
  flow‑cell. Some samples required re-runs and were pooled before
  alignment/variant calling to maximise depth.
- This initially requires UMIs to be parsed and extracted into read-names via 
`process_fastq_umi_for_merging.sh`
- Samples can then be concatenated together and returned to Austin Path for DRAGEN processing i.e.  
`cat Sample1_Batch1_R1.fastq.gz Sample1_Batch2_R1_fastq.gz > Sample1_Batch1_Batch2_R1.fastq.gz`  
`cat Sample1_Batch1_R2.fastq.gz Sample1_Batch2_R2_fastq.gz > Sample1_Batch1_Batch2_R2.fastq.gz`

---

### Duplex sequencing explanation

Further information:

<https://help.dragen.illumina.com/product-guide/dragen-v4.4/dragen-dna-pipeline/unique-molecular-identifiers>

Illumina tech support explanation of XV/XW:

> "A read pair of a duplex family will have 2 scenarios of simplex
> (R1 forward and R2 reverse, or R1 reverse and R2 forward).  One
> strand will be randomly selected for the consensus of the collapsed
> read.  The XV represents the number of raw read pairs matching to the
> selected consensus strand and XW provides the number of raw read
> pairs opposite to consensus strand.  In some cases, the opposite
> strand may have larger number of the read supports, so XW can be
> larger than XV.  If consensus is simplex, then there's no duplex
> pair, so XW will be 0."

---

## Data resources

- `/vast/projects/reidj-project/cosmic/Cosmic_GenomeScreensMutant_v102_GRCh38.vcf.gz`
- `/vast/projects/reidj-project/clinvar/clinvar_20250330.vcf.gz`
- Singularity/Apptainer containers:
  - `/vast/projects/reidj-project/containers/cadd-scoring_latest.sif`
  - `/stornext/Bioinf/data/lab_bahlo/users/reid.j/analysis/Apptainer/spliceai.img`
- External annotations:
  - `/stornext/Bioinf/data/lab_bahlo/users/reid.j/cadd/CADD-scripts-master/data/annotations`
  - `/stornext/Bioinf/data/lab_bahlo/users/reid.j/spliceai/gencode.v43.canonical.annotation.txt`

## TODO

- [ ] 

---
