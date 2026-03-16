#!/usr/bin/env Rscript
# process_variant.R
# This script processes one VCF file and one corresponding variant data TSV file.
# It reads the input files from the command line and writes a processed output as an RDS file.

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(stringr)
  library(dplyr)
  library(tidyr)
})

# Check and collect command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript filter_austin_panel_variants.R <vcf> <vaf_csv> <output_rds> <project_dir> \
       <min_duplex> <min_simplex> <min_duplex_with_simplex> <min_simplex_with_duplex>")
}

vcf_file                <- args[1]
variant_file            <- args[2]
output_rds              <- args[3]
project_dir             <- args[4]
min_duplex              <- as.integer(args[5])
min_simplex             <- as.integer(args[6])
min_duplex_with_simplex <- as.integer(args[7])
min_simplex_with_duplex <- as.integer(args[8])


# Use project_dir to build paths to gene lists
ras <- read.table(file.path(project_dir, "pipeline_files/gene_lists/20240923_ras_pathway_genes_namesonly.tsv"),
                  sep = "\t", header = FALSE)
mtor <- read.table(file.path(project_dir, "pipeline_files/gene_lists/20240923_mtor_pathway_genes_namesonly.tsv"),
                   sep = "\t", header = FALSE)
mcd_panelapp <- read.table(file.path(project_dir, "pipeline_files/gene_lists/MCD_Superpanel_PanelApp_GeneNames_v4.65.tsv"),
                            sep = "\t", header = FALSE)
g4e <- read.table(file.path(project_dir, "pipeline_files/gene_lists/EpilepsyGenes_v2024-09.tsv"),
                  sep = "\t", header = TRUE)
colnames(g4e)[7] <- 'Phenotype'

gene_info <- as.data.frame(c(ras$V1, mtor$V1, mcd_panelapp$V1, g4e$Gene))
colnames(gene_info) <- 'Gene'
gene_info <- unique(gene_info)

gene_info$Pathway <- ''

for (k in 1:nrow(gene_info)){
  gene <- gene_info[k,'Gene']
  pathways <- character(0)  # Initialize an empty character vector to store pathways
  
  if (gene %in% ras$V1){
    pathways <- c(pathways, 'RAS')
  }
  if (gene %in% mtor$V1){
    pathways <- c(pathways, 'MTOR')
  }
  if (gene %in% g4e$Gene){
    pathways <- c(pathways, 'G4E')
  }
  if (gene %in% mcd_panelapp$V1){
    pathways <- c(pathways, 'MCD_PanelApp')
  }
  
  # Concatenate the pathways if multiple conditions are true
  gene_info[k,'Pathway'] <- paste(pathways, collapse = ', ')
}

# Set reference genome
ref <- "hg38"

# Helper function to reformat variants for Franklin URL
reformat_variants <- function(variants) {
  base_url <- "https://franklin.genoox.com/clinical-db/variant/snp/"
  urls <- paste(base_url, gsub(":", "-", variants), "-hg38", sep = "")
  return(urls)
}

# Check if VCF has variants BEFORE reading
n_variants <- as.integer(system(
    sprintf("zgrep -v '^#' %s | wc -l", vcf_file),
    intern = TRUE
))

if (n_variants == 0) {
  cat("VCF file is empty. Saving empty data frame.\n")
  
  # Initialize empty data frame with correct structure
  vcf_info <- data.frame(
    CHR = character(0),
    POS = numeric(0),
    REF = character(0),
    ALT = character(0),
    Variant = character(0),
    Consequence = character(0),
    IMPACT = character(0),
    SYMBOL = character(0),
    Gene = character(0),
    Feature = character(0),
    BIOTYPE = character(0),
    EXON = character(0),
    INTRON = character(0),
    HGVSc = character(0),
    HGVSp = character(0),
    HGVSg = character(0),
    Existing_variation = character(0),
    DISTANCE = character(0),
    FLAGS = character(0),
    VARIANT_CLASS = character(0),
    CADD_Score = numeric(0),
    SIFT = character(0),
    PolyPhen = character(0),
    CLIN_SIG = character(0),
    AC_gnomad_exomes_4.0 = numeric(0),
    AC_gnomad_genomes_4.0 = numeric(0),
    AC_gnomad_total_4.0 = numeric(0),
    nhomalt_gnomad_exomes_4.0 = numeric(0),
    nhomalt_gnomad_genomes_4.0 = numeric(0),
    nhomalt_gnomad_total_4.0 = numeric(0),
    COSMIC_Sample_Count = numeric(0),
    REF_AD = numeric(0),
    ALT_AD = numeric(0),
    VAF = numeric(0),
    BAM_REF = numeric(0),
    BAM_ALT = numeric(0),
    BAM_OTHER = numeric(0),
    DUPLEX_REF = numeric(0),
    DUPLEX_ALT = numeric(0),
    DUPLEX_OTHER = numeric(0),
    SIMPLEX_REF = numeric(0),
    SIMPLEX_ALT = numeric(0),
    SIMPLEX_OTHER = numeric(0),
    Gene.1 = character(0),
    Pathway = character(0),
    Inheritance = character(0),
    Phenotype = character(0),
    Franklin_URL = character(0),
    stringsAsFactors = FALSE
  )
  
  saveRDS(vcf_info, file = output_rds)
  cat("Empty data saved to", output_rds, "\n")
  quit(save = "no", status = 0)
}

# Now safe to read VCF - we know it has variants
vcf <- readVcf(vcf_file, ref)
sample_names <- colnames(geno(vcf)$GT)
if(length(sample_names) == 0){
  stop("No sample found in VCF file.")
}
sample_name <- sample_names[1]
  
# Create a basic data frame from the VCF
vcf_info <- data.frame(
  CHR = as.character(seqnames(vcf)),
  POS = start(vcf),
  REF = as.character(ref(vcf)),
  ALT = sapply(alt(vcf), function(x) as.character(x[1])),
  stringsAsFactors = FALSE
)

if (nrow(vcf_info) == 0) {
  # If VCF is empty, return an empty df with the correct column names and save as RDS
  vcf_info <- empty_vcf_info
  saveRDS(vcf_info, file = output_rds)
  cat("VCF is empty. Data saved to", output_rds, "\n")
} else {
  vcf_info$Variant <- paste(vcf_info$CHR, vcf_info$POS, vcf_info$REF, vcf_info$ALT, sep = ":")

# Process AD field for allele depths if present
ad <- geno(vcf)$AD
if (!is.null(ad)) {
  ref_alt_counts <- ad[, 1]
  ref_ad_values <- sapply(ref_alt_counts, function(counts) if (is.na(counts[1])) NA else counts[1])
  alt_ad_values <- sapply(ref_alt_counts, function(counts) if (is.na(counts[2])) NA else counts[2])
  
  if (nrow(vcf_info) == 0) {
    vcf_info$REF_AD <- numeric(0)
    vcf_info$ALT_AD <- numeric(0)
    vcf_info$VAF <- numeric(0)
  } else {
    vcf_info$REF_AD <- ref_ad_values
    vcf_info$ALT_AD <- alt_ad_values
    vcf_info$VAF <- alt_ad_values / (ref_ad_values + alt_ad_values)
  }
}

# Extract the INFO fields
info_fields <- data.frame(info(vcf))
vcf_info <- cbind(vcf_info, info_fields)

# Read in the variant data file generated by the Python script
variant_data <- read.table(variant_file, header = TRUE, sep = ",", stringsAsFactors = FALSE,
                           colClasses = c("REF" = "character", "ALT" = "character"))

# Ensure required columns exist even when vcf_info is empty
if (nrow(vcf_info) == 0) {
  vcf_info <- data.frame(
    CHR = character(0),
    POS = integer(0),
    REF = character(0),
    ALT = character(0),
    Variant = character(0),
    VAF = numeric(0),
    REF_AD = numeric(0),
    ALT_AD = numeric(0),
    stringsAsFactors = FALSE
  )
}

# Merge the VCF info with the variant data using common columns
vcf_info <- merge(vcf_info, variant_data, by = c("CHR", "POS", "REF", "ALT"), all.x = TRUE)

#Replicate variant rows have multiple VEP entries (i.e. variants that impact multiple genes)
vcf_info <- vcf_info %>% 
  unnest(CSQ) %>% 
  distinct(CSQ, .keep_all = TRUE)

##VEP CSQ##
#VEP consequences are split into individual fields
if ("CSQ" %in% colnames(vcf_info)){
  print('Processing Ensembl VEP CSQ')
  csq = info(header(vcf))['CSQ',]$Description
  csq_col = unlist(strsplit(csq, ':'))[2]
  csq_fields = unlist(strsplit(csq_col, '\\|'))
  csq_data <- data.frame(matrix(nrow = nrow(vcf_info), ncol = length(csq_fields)))
  colnames(csq_data) <- csq_fields
  for (i in 1:nrow(vcf_info)){
    transcript <- unlist(vcf_info[i,]$CSQ)
    fields <- unlist(strsplit(transcript, split = '\\|'))
    for (col in 1:length(csq_fields)){
      csq_data[i,col] <- fields[col]}
  }
  vcf_info <- cbind(csq_data, vcf_info)}

## gnomAD ##
if ("AC_gnomad_exomes_4.0" %in% colnames(vcf_info)){
  print('Processing AC_gnomad_exomes_4.0')
  for (i in 1:length(vcf_info$AC_gnomad_exomes_4.0)){
    entry = vcf_info$AC_gnomad_exomes_4.0[i]
    if (length(unlist(entry)) > 1) {
      max = 0
      for (j in 1:length(unlist(entry))){
        if (!is.na(unlist(entry)[j])){
          count = as.numeric(unlist(entry)[j])
          if (count > max){
            max = count}
        }
        entry = max}
    } else if (is.na(entry)){
      entry = 0
    }
    vcf_info$AC_gnomad_exomes_4.0[i] = as.numeric(entry)
  }
  vcf_info$AC_gnomad_exomes_4.0 <- as.numeric(vcf_info$AC_gnomad_exomes_4.0)
}

if ("AC_gnomad_genomes_4.0" %in% colnames(vcf_info)){
  print('Processing AC_gnomad_genomes_4.0')
  for (i in 1:length(vcf_info$AC_gnomad_genomes_4.0)){
    entry = vcf_info$AC_gnomad_genomes_4.0[i]
    if (length(unlist(entry)) > 1) {
      max = 0
      for (j in 1:length(unlist(entry))){
        if (!is.na(unlist(entry)[j])){
          count = as.numeric(unlist(entry)[j])
          if (count > max){
            max = count}
        }
        entry = max}
    } else if (is.na(entry)){
      entry = 0
    }
    vcf_info$AC_gnomad_genomes_4.0[i] = as.numeric(entry)
  }
  vcf_info$AC_gnomad_genomes_4.0 <- as.numeric(vcf_info$AC_gnomad_genomes_4.0)
}

vcf_info$AC_gnomad_total_4.0 <- as.numeric(vcf_info$AC_gnomad_exomes_4.0) + as.numeric(vcf_info$AC_gnomad_genomes_4.0)

if ("nhomalt_gnomad_exomes_4.0" %in% colnames(vcf_info)){
  print('Processing nhomalt_gnomad_exomes_4.0')
  for (i in 1:length(vcf_info$nhomalt_gnomad_exomes_4.0)){
    entry = vcf_info$nhomalt_gnomad_exomes_4.0[i]
    if (length(unlist(entry)) > 1) {
      max = 0
      for (j in 1:length(unlist(entry))){
        if (!is.na(unlist(entry)[j])){
          count = as.numeric(unlist(entry)[j])
          if (count > max){
            max = count}
        }
        entry = max}
    } else if (is.na(entry)){
      entry = 0
    }
    vcf_info$nhomalt_gnomad_exomes_4.0[i] = as.numeric(entry)
  }
  vcf_info$nhomalt_gnomad_exomes_4.0 <- as.numeric(vcf_info$nhomalt_gnomad_exomes_4.0)
}

if ("nhomalt_gnomad_genomes_4.0" %in% colnames(vcf_info)){
  print('Processing nhomalt_gnomad_genomes_4.0')
  for (i in 1:length(vcf_info$nhomalt_gnomad_genomes_4.0)){
    entry = vcf_info$nhomalt_gnomad_genomes_4.0[i]
    if (length(unlist(entry)) > 1) {
      max = 0
      for (j in 1:length(unlist(entry))){
        if (!is.na(unlist(entry)[j])){
          count = as.numeric(unlist(entry)[j])
          if (count > max){
            max = count}
        }
        entry = max}
    } else if (is.na(entry)){
      entry = 0
    }
    vcf_info$nhomalt_gnomad_genomes_4.0[i] = as.numeric(entry)
  }
  vcf_info$nhomalt_gnomad_genomes_4.0 <- as.numeric(vcf_info$nhomalt_gnomad_genomes_4.0)
}

vcf_info$nhomalt_gnomad_total_4.0 <- as.numeric(vcf_info$nhomalt_gnomad_exomes_4.0) + as.numeric(vcf_info$nhomalt_gnomad_genomes_4.0)

#COSMIC
vcf_info$COSMIC_Sample_Count <- as.numeric(vcf_info$COSMIC_Sample_Count)

# Convert all SpliceAI columns to numeric
vcf_info$SpliceAI_DS_AG <- as.numeric(vcf_info$SpliceAI_DS_AG)
vcf_info$SpliceAI_DS_AL <- as.numeric(vcf_info$SpliceAI_DS_AL)
vcf_info$SpliceAI_DS_DG <- as.numeric(vcf_info$SpliceAI_DS_DG)
vcf_info$SpliceAI_DS_DL <- as.numeric(vcf_info$SpliceAI_DS_DL)
vcf_info$SpliceAI_Max <- as.numeric(vcf_info$SpliceAI_Max)

#AlphaMissense
colnames(vcf_info)[colnames(vcf_info) == "am_class"] <- "AlphaMissense_Class"
colnames(vcf_info)[colnames(vcf_info) == "am_pathogenicity"] <- "AlphaMissense_Score"

# Match VCF genes with Epilepsy Gene List genes
matches <- match(vcf_info$SYMBOL, gene_info$Gene)

# Subset gene_info based on the matching rows
gene_summary <- gene_info[matches,]

# Add the matched gene summary to vcf_info
vcf_info <- cbind(vcf_info, gene_summary)

# Merge with g4e to add Inheritance and Phenotype without duplicating columns
vcf_info <- merge(
  vcf_info, 
  g4e[, c("Gene", "Inheritance", "Phenotype")],  # Include "Gene" here for the join
  by.x = "SYMBOL", 
  by.y = "Gene", 
  all.x = TRUE, 
  suffixes = c("", "_g4e")
)

#Remove repeat entries (i.e. same variant across multiple transcripts). 
# Find all duplicated variants and their corresponding rows
rownames(vcf_info) <- seq(nrow(vcf_info))

duplicate_variant_index_1 <- which(duplicated(vcf_info$Variant) | duplicated(vcf_info$Variant, fromLast = TRUE))

#First, remove any duplicates that have a blank SYMBOL column (i.e. no gene name)
vcf_info <- subset(vcf_info, !((rownames(vcf_info) %in% duplicate_variant_index_1) & (SYMBOL == "")))

# Reformat variant URLs for Franklin lookup
vcf_info$Franklin_URL <- reformat_variants(vcf_info$Variant)

vcf_info <- vcf_info %>%
  dplyr::select(
    CHR, POS, REF, ALT, Variant, Consequence, IMPACT, SYMBOL, Gene, Feature, 
    BIOTYPE, EXON, INTRON, HGVSc, HGVSp, HGVSg, Existing_variation, DISTANCE, 
    FLAGS, VARIANT_CLASS, CADD_Score, SIFT, AlphaMissense_Class, AlphaMissense_Score, REVEL,
    PolyPhen, SpliceAI_DS_AG, SpliceAI_DS_AL, SpliceAI_DS_DG, SpliceAI_DS_DL, SpliceAI_Max,
    ClinVarSIG, ClinVarSIGCONF, ClinVarGene, ClinVarDN, 
    AC_gnomad_exomes_4.0, AC_gnomad_genomes_4.0, AC_gnomad_total_4.0, 
    nhomalt_gnomad_exomes_4.0, nhomalt_gnomad_genomes_4.0, nhomalt_gnomad_total_4.0, COSMIC_Sample_Count, 
    REF_AD, ALT_AD, VAF, BAM_REF, BAM_ALT, BAM_OTHER, DUPLEX_REF, DUPLEX_ALT, DUPLEX_OTHER, 
    SIMPLEX_REF, SIMPLEX_ALT, SIMPLEX_OTHER, Gene.1, Pathway, Inheritance, Phenotype, Franklin_URL
  )

# Filter vcf_info based on duplex/simplex coverage
vcf_info <- vcf_info %>%
  filter(
    DUPLEX_ALT >= min_duplex |
    SIMPLEX_ALT >= min_simplex |
    (DUPLEX_ALT >= min_duplex_with_simplex & SIMPLEX_ALT >= min_simplex_with_duplex)
  )

# Save the processed data frame as an RDS file
saveRDS(vcf_info, file = output_rds)
cat("Processing complete. Data saved to", output_rds, "\n")
}