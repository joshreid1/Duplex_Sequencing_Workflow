#!/usr/bin/env Rscript
library(vcfR)           # vcfR for reading and manipulating VCF files
library(ggplot2)
library(reshape2)
library(dndscv)
library(dplyr)

?dndscv

gene_panel <- scan(text = "
AKT1 ALK MAP2K1 ATP2A1 ACVR2A CTNNB1 IDH2 NCOA4 SDC4
AKT3 BARD1 MAP3K3 GNA11 AIM2 DAXX IL6ST SDHA
DEPDC5 BRAF MDM2 GNAQ APC DEPDC5 KBTBD4 SDHB
MTOR BRCA1 MET GRIN2C ARID1A DICER1 KDM6A PALB2 SDHC
NPRL2 BRCA2 MYB KLHL22 ASTE1 DROSHA KEAP1 PBRM1 SDHD
NPRL3 BRIP1 SF3B1 ATM EIF1AX KIF5B PDGFRA SLC34A2
PIK3CA CBL NF1 SLC35A2 ATR ELOC KIT PMS2 SMAD4
PIK3R1 CCND1 NF2 ATRX EML4 KMT2C POLE SMARCA4
PTEN CDK4 NRAS BAP1 EPCAM KMT2D PPFIA4 TAF1B
RHEB EGFR PRKCA EZR MARCKS PRKAR1A TEK
TSC1 ERBB2 PTPN11 FANCL MEN1 PRKD1 TERT
TSC2 FANCA RAF1 CCDC6 FH MITF PTHLH TFEB
FGFR1 RB1 CD74 FOXL2 MLH1 TGFBR2 FGFR2
RIT1 CDK12 GNAS RAD51B FGFR3 ROS1 CDKN2A
H3-3A MRE11 RAD51C HRAS SOS1 CDKN2B H3C2 MSH2
RAD54L KRAS STK11 CHEK1 HNF1A MSH6 RET LZTR1
TP53 CHEK2 IDH1 NBN RNF43
", what = "")

#The following input gene names are not in the RefCDS database: NR21, NR24, MYC-N, BAT25, BAT26, RAD51A, MONO27
length(gene_panel)

### Test Data ###

data("dataset_simbreast", package="dndscv")

####################
setwd("~/vast_scratch/austin_panel/dndscv/")
sample_data <- read.table("/vast/projects/reidj-project/Austin_Panel/austin_panel_id_groups_plus_vcfpath.tsv", sep ="\t", header = TRUE)

sample_data[1,]

#> sample_data[1,]
#Group Sample_ID  Batch
#1  MTLE   38958-1 DOUBLE
#Filepath
#1 /vast/scratch/users/reid.j/austin_panel/Oct_20/2025-10-15_Merge_Dragen/38958-1-EST007-EST026_ds.e436ac048fcd44a8a5680c5680a4a650/38958-1-EST007-EST026.hard-filtered.vcf.gz

#test("~/vast_scratch/austin_panel/EST032_Dragen/50758-2_ds.00a1dac718b843d2be266f5d4d6d2590/50758-2.sv.vcf.gz"

mtle_orignal <- sample_data[sample_data$Group == "MTLE", 'Filepath']
mtle <- mtle_orignal[!grepl("EST035|EST036", mtle_orignal)]
mnd <- sample_data[sample_data$Group == "MND", 'Filepath']

mnd

#Test of similar cohort size
#mtle <- mtle[41:81]

##Test on random assignment
#mtle1 <- mtle[seq(1, length(mtle), by = 2)]
#mtle2 <- mtle[seq(2, length(mtle), by = 2)]
#mtle <- mtle1
#mnd <- mtle2


#                         min_total_depth = 50, 
#                        max_total_depth = NA, 
#                         min_alt_depth_snp = 3,
#                         min_alt_depth_indel = 3,
#                         max_vaf = 0.25, 
#                         min_vaf = 0.002) {


# Helper: Filter variants by allelic depth (FORMAT/AD) with separate thresholds for SNPs/indels
filter_by_ad <- function(vcf_object, 
                         min_total_depth = NA, 
                         max_total_depth = NA, 
                         min_alt_depth_snp = 3,
                         min_alt_depth_indel = 3,
                         max_vaf = 0.25, 
                         min_vaf = NA) {
  
  gt <- extract.gt(vcf_object, element = "AD") # matrix, samples as columns
  
  # AD in GATK-like VCFs: comma-separated, e.g., "8,2"
  ad_split <- strsplit(gt[,1], ",")   # assumes one sample per VCF
  total_depth <- sapply(ad_split, function(x) sum(as.numeric(x)))
  alt_depth  <- sapply(ad_split, function(x) as.numeric(x[2]))
  vaf        <- alt_depth / total_depth
  
  # Determine variant type (SNP vs indel)
  ref_alleles <- as.character(vcf_object@fix[, "REF"])
  alt_alleles <- as.character(vcf_object@fix[, "ALT"])
  is_snp <- nchar(ref_alleles) == 1 & nchar(alt_alleles) == 1
  
  # Initialize keep vector - start with TRUE for all variants
  keep <- !is.na(alt_depth) & !is.na(total_depth)
  
  # Apply filters only if parameter is not NA
  if (!is.na(min_total_depth)) {
    keep <- keep & total_depth >= min_total_depth
  }
  
  if (!is.na(max_total_depth)) {
    keep <- keep & total_depth <= max_total_depth
  }
  
  # Apply variant-type-specific min_alt_depth thresholds
  if (!is.na(min_alt_depth_snp)) {
    keep <- keep & ((!is_snp) | (alt_depth >= min_alt_depth_snp))
  }
  
  if (!is.na(min_alt_depth_indel)) {
    keep <- keep & (is_snp | (alt_depth >= min_alt_depth_indel))
  }
  
  if (!is.na(max_vaf)) {
    keep <- keep & vaf <= max_vaf
  }
  
  if (!is.na(min_vaf)) {
    keep <- keep & vaf >= min_vaf
  }
  
  return(keep)
}


# Function to remove variants with multiple alternative alleles
remove_mult_alts = function(vcf_df) {
  mult_alts = grepl(",", vcf_df$ALT)
  nr_mult_alts = sum(mult_alts)
  if (nr_mult_alts > 0){
    vcf_df = vcf_df[!mult_alts, ]
    warning(paste0("Removed ", nr_mult_alts, " variant(s), because they had multiple alternative alleles."), call. = F)
  }
  return(vcf_df)
}

# Modified: Plot VAF distribution from vcfR object(s) or list of vcfR objects
plot_vaf_distribution <- function(vcf_object, sample_col = 1, 
                                  title = "VAF Distribution",
                                  binwidth = 0.01) {
  require(ggplot2)
  require(vcfR)
  
  # Handle list of vcfR objects (for cohort-level plotting)
  if (is.list(vcf_object) && !inherits(vcf_object, "vcfR")) {
    # Combine multiple vcfR objects
    # Filter out NULL entries (from samples that didn't pass filters)
    vcf_object <- vcf_object[!sapply(vcf_object, is.null)]
    
    if (length(vcf_object) == 0) {
      warning("No valid vcfR objects to plot")
      return(NULL)
    }
    
    # Combine all vcfR objects using rbind
    vcf_object <- do.call(rbind, vcf_object)
  }
  
  # Extract AD field (allelic depth)
  gt <- extract.gt(vcf_object, element = "AD")
  
  # Split comma-separated AD values (e.g., "8,2")
  ad_split <- strsplit(gt[, sample_col], ",")
  
  # Calculate total depth, alt depth, and VAF
  total_depth <- sapply(ad_split, function(x) sum(as.numeric(x)))
  alt_depth   <- sapply(ad_split, function(x) as.numeric(x[2]))
  vaf         <- alt_depth / total_depth
  
  # Create data frame, removing NAs
  vaf_df <- data.frame(VAF = vaf[!is.na(vaf)])
  
  # Create histogram with density overlay
  p <- ggplot(vaf_df, aes(x = VAF)) +
    geom_histogram(aes(y = after_stat(density)), 
                   binwidth = binwidth, 
                   fill = "steelblue", 
                   color = "black", 
                   alpha = 0.7) +
    geom_density(color = "red", linewidth = 1) +
    labs(title = title,
         x = "Variant Allele Frequency (VAF)",
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# Modified: Load mutations with cohort-level VAF plotting
load_mutations <- function(sampleslist, group_name) {
  mutations <- data.frame()
  
  # Initialize lists to accumulate vcfR objects across the cohort
  all_unfiltered_vcf <- list()
  all_filtered_vcf <- list()
  
  for(i in 1:length(sampleslist)){
    unfiltered_vcf_object <- read.vcfR(sampleslist[i], verbose = FALSE)
    
    # Store unfiltered vcfR object
    all_unfiltered_vcf[[i]] <- unfiltered_vcf_object
    
    # AD pre-filter (before getFIX, extract.info, etc)
    ad_keep <- filter_by_ad(unfiltered_vcf_object)
    if (!any(ad_keep)) {
      all_filtered_vcf[[i]] <- NULL  # Mark as NULL if nothing passes
      next # skip if nothing passes
    }
    
    # get only rows passing AD filter for subsequent processing
    vcf_object <- unfiltered_vcf_object[ad_keep,]
    
    # Store filtered vcfR object
    all_filtered_vcf[[i]] <- vcf_object
    
    # Extract sample ID from filename (adjust if needed)
    SampleID <- basename(sampleslist[i])
    print(paste0("Processing ", group_name, " sample: ", SampleID))
    
    vcf_fix <- getFIX(vcf_object, getINFO = FALSE)
    vcf_fix <- as.data.frame(vcf_fix, stringsAsFactors = FALSE)
    
    vcf_info <- extract.info(vcf_object, element = NULL)
    vcf_df <- cbind(vcf_fix, vcf_info)
    
    vcf_df <- remove_mult_alts(vcf_df)
    vcf_df <- vcf_df[vcf_df$FILTER == "PASS", ]
    
    chrom <- as.character(sub("chr", "", vcf_df$CHROM))
    
    mutations_sample <- data.frame(
      sampleID = SampleID,
      chr = chrom,
      pos = as.numeric(vcf_df$POS),
      ref = vcf_df$REF,
      mut = vcf_df$ALT,
      stringsAsFactors = FALSE
    )
    
    # Remove rows with any NA values in key columns
    mutations_sample <- mutations_sample %>%
      filter(!is.na(chr), !is.na(pos), !is.na(ref), !is.na(mut))
    
    # Optional: Print warning if NAs were found and removed
    n_removed <- nrow(mutations_sample) - sum(complete.cases(mutations_sample))
    if (n_removed > 0) {
      print(paste0("  Warning: Removed ", n_removed, " variants with NA values from ", SampleID))
    }
    
    mutations <- rbind(mutations, mutations_sample)
  }
  
  # Plot cohort-level VAF distributions after processing all samples
  print(paste0("Plotting cohort-level VAF distributions for ", group_name))
  
  # Count total variants for titles
  n_unfiltered <- sum(sapply(all_unfiltered_vcf, function(x) if(!is.null(x)) nrow(x) else 0))
  n_filtered <- sum(sapply(all_filtered_vcf, function(x) if(!is.null(x)) nrow(x) else 0))
  
  p_unfiltered <- plot_vaf_distribution(
    all_unfiltered_vcf, 
    title = paste0(group_name, " - Unfiltered (n=", n_unfiltered, " variants)")
  )
  
  p_filtered <- plot_vaf_distribution(
    all_filtered_vcf, 
    title = paste0(group_name, " - Filtered (n=", n_filtered, " variants)")
  )
  
  # Print plots explicitly to display them
  print(p_unfiltered)
  print(p_filtered)
  
  return(mutations)
}



####################
# Analyze both cohorts
####################

cat("\n=== Loading MTLE cohort ===\n")
mutations_mtle <- load_mutations(mtle, "MTLE")
cat(paste0("\nMTLE: ", length(unique(mutations_mtle$sampleID)), " samples, ", nrow(mutations_mtle), " variants\n"))
head(mutations_mtle)
sample_counts_mtle <- table(mutations_mtle$sampleID)
summary(as.numeric(sample_counts_mtle))

cat("\n=== Loading MND cohort ===\n")
mutations_mnd <- load_mutations(mnd, "MND")
cat(paste0("\nMND: ", length(unique(mutations_mnd$sampleID)), " samples, ", nrow(mutations_mnd), " variants\n"))
sample_counts_mnd <- table(mutations_mnd$sampleID)
summary(as.numeric(sample_counts_mnd))

# Create a combined data frame for plotting
counts_df <- data.frame(
  sample_id = c(names(sample_counts_mtle), names(sample_counts_mnd)),
  variant_count = c(as.numeric(sample_counts_mtle), as.numeric(sample_counts_mnd)),
  cohort = c(rep("MTLE", length(sample_counts_mtle)), 
             rep("MND", length(sample_counts_mnd)))
)

#Boxplot of filtered variants per sample
ggplot(counts_df, aes(x = cohort, y = variant_count, fill = cohort)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  labs(x = "Cohort", y = "Variant count per sample",
       title = "Variant count distribution by cohort") +
  theme_bw()

head(counts_df)

# Load the age data
age_data <- read.table("/stornext/General/scratch/GP_Transfer/reid.j/austin_panel_age_data.tsv", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Clean sample_id in counts_df: remove .hard-filtered.vcf.gz suffix
counts_df$sample_id_clean <- gsub("\\.hard-filtered\\.vcf\\.gz$", "", counts_df$sample_id)

# Remove -* suffix from both dataframes for matching
counts_df$sample_id_base <- gsub("-.*$", "", counts_df$sample_id_clean)
age_data$sample_id_base <- gsub("-.*$", "", age_data$sample_id)

# Merge on the base sample_id
merged_df <- merge(counts_df, age_data, by = "sample_id_base", suffixes = c("", "_age"))

# Display the merged data
cat("\nMerged data summary:\n")
print(head(merged_df))
cat(paste0("\nTotal samples with age data: ", nrow(merged_df), "\n"))

ggplot(merged_df, aes(x = age, y = variant_count, color = cohort)) +
  geom_point(size = 3, alpha = 0.7) +
  #scale_y_log10() +
  labs(x = "Age (years)", 
       y = "Somatic variant count",
       title = "Variant count vs. Age by cohort",
       color = "Cohort",
       fill = "Cohort") +
  theme_bw() +
  theme(legend.position = "right")

head(merged_df)
head(depth_data)

# Calculate median depth per sample
depth_summary <- depth_data %>%
  group_by(sampleID) %>%
  summarise(median_depth = median(total_depth, na.rm = TRUE))

# Extract sample_id_base to match merged_df structure
depth_summary <- depth_summary %>%
  mutate(sample_id_clean = gsub("\\.hard-filtered\\.vcf\\.gz$", "", sampleID),
         sample_id_base = sub("-.*", "", sample_id_clean))

# Merge with existing merged_df
merged_df_plus <- merge(merged_df, 
                   depth_summary[, c("sample_id_base", "median_depth")], 
                   by = "sample_id_base", 
                   all.x = TRUE)

# Primary normalization: variants per 100x coverage
merged_df_plus$variants_per_100x <- (merged_df_plus$variant_count / merged_df_plus$median_depth) * 100

library(patchwork)

# Compare raw vs. normalized
pA <- ggplot(merged_df_plus, aes(x = age, y = variant_count, color = cohort)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  labs(x = "Age (years)", y = "Variant count", 
       title = "Raw variant count vs. Age") +
  theme_bw()

pB <- ggplot(merged_df_plus, aes(x = age, y = variants_per_100x, color = cohort)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  labs(x = "Age (years)", y = "Variants per 100x coverage",
       title = "Depth-normalized variants vs. Age") +
  theme_bw()

# QC: Check if depth varies with age
pC<- ggplot(merged_df_plus, aes(x = age, y = median_depth, color = cohort)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "Sequencing depth QC") +
  theme_bw()

pA / pB / pC

# Optional: Print correlation statistics by cohort
cat("\n=== Correlation statistics ===\n")
for (cohort_name in unique(merged_df$cohort)) {
  cohort_data <- merged_df[merged_df$cohort == cohort_name, ]
  if (nrow(cohort_data) > 2) {
    cor_test <- cor.test(cohort_data$age, cohort_data$variant_count)
    cat(paste0(cohort_name, ": r = ", round(cor_test$estimate, 3), 
               ", p = ", format.pval(cor_test$p.value, digits = 3), 
               " (n = ", nrow(cohort_data), ")\n"))
  }
}

####################
# Run dNdScv on both cohorts
####################

cat("\n=== Running dNdScv on MTLE cohort ===\n")
dnds_mtle <- dndscv(mutations_mtle, refdb="hg38", 
                    gene_list = gene_panel,
                    max_muts_per_gene_per_sample = Inf,
                    max_coding_muts_per_sample = Inf,
                    constrain_wnon_wspl = TRUE,
                    outmats = TRUE,
                    mingenecovs = 50)

dnds_mtle <- dndscv(mutations_mtle, refdb="hg38", 
                    gene_list = gene_panel,
                    max_muts_per_gene_per_sample = Inf,
                    max_coding_muts_per_sample = Inf,
                    outmats = TRUE)

cat("\n=== Running dNdScv on MND cohort ===\n")
dnds_mnd <- dndscv(mutations_mnd, refdb="hg38",
                   gene_list = gene_panel,
                   max_muts_per_gene_per_sample = Inf,
                   max_coding_muts_per_sample = Inf,
                   constrain_wnon_wspl = TRUE,
                   outmats = TRUE)

dnds_mnd <- dndscv(mutations_mnd, refdb="hg38", 
                    gene_list = gene_panel,
                    max_muts_per_gene_per_sample = Inf,
                    max_coding_muts_per_sample = Inf,
                    outmats = TRUE)

####################
# Extract and compare results
####################

sel_mtle <- dnds_mtle$sel_cv
sel_mnd <- dnds_mnd$sel_cv

# Add group identifier
sel_mtle$group <- "MTLE"
sel_mnd$group <- "MND"

# Global dN/dS ratios
global_mtle <- dnds_mtle$globaldnds
global_mnd <- dnds_mnd$globaldnds

cat("\n=== Global dN/dS ratios ===\n")
cat("\nMTLE:\n")
print(global_mtle)
cat("\nMND:\n")
print(global_mnd)

# Combine global results for comparison
global_comparison <- data.frame(
  group = c("MTLE", "MND"),
  wmis = c(global_mtle$mle[global_mtle$name == "wmis"], 
           global_mnd$mle[global_mnd$name == "wmis"]),
  wnon = c(global_mtle$mle[global_mtle$name == "wnon"], 
           global_mnd$mle[global_mnd$name == "wnon"]),
  wspl = c(global_mtle$mle[global_mtle$name == "wspl"], 
           global_mnd$mle[global_mnd$name == "wspl"])
)

write.csv(global_comparison, "global_dnds_comparison.csv", row.names = FALSE)

####################
# Statistical comparison of genes
####################

# Merge gene-level results
gene_comparison <- merge(
  sel_mtle[, c("gene_name", "wmis_cv", "wnon_cv", "wspl_cv", "qglobal_cv")],
  sel_mnd[, c("gene_name", "wmis_cv", "wnon_cv", "wspl_cv", "qglobal_cv")],
  by = "gene_name",
  suffixes = c("_MTLE", "_MND")
)

# Calculate differences
gene_comparison$delta_wmis <- gene_comparison$wmis_cv_MTLE - gene_comparison$wmis_cv_MND
gene_comparison$delta_wnon <- gene_comparison$wnon_cv_MTLE - gene_comparison$wnon_cv_MND
gene_comparison$delta_wspl <- gene_comparison$wspl_cv_MTLE - gene_comparison$wspl_cv_MND

# Identify genes significant in either cohort
gene_comparison$sig_MTLE <- gene_comparison$qglobal_cv_MTLE < 0.05
gene_comparison$sig_MND <- gene_comparison$qglobal_cv_MND < 0.05
gene_comparison$sig_either <- gene_comparison$sig_MTLE | gene_comparison$sig_MND

write.csv(gene_comparison, "gene_level_comparison.csv", row.names = FALSE)

gene_comparison[gene_comparison$sig_MTLE == TRUE, ]

# Summary of significant genes
cat("\n=== Significant genes (q < 0.1) ===\n")
cat(paste0("\nMTLE only: ", sum(gene_comparison$sig_MTLE & !gene_comparison$sig_MND), " genes\n"))
cat(paste0("MND only: ", sum(gene_comparison$sig_MND & !gene_comparison$sig_MTLE), " genes\n"))
cat(paste0("Both cohorts: ", sum(gene_comparison$sig_MTLE & gene_comparison$sig_MND), " genes\n"))

if(sum(gene_comparison$sig_either) > 0) {
  sig_genes_table <- gene_comparison[gene_comparison$sig_either, 
                                     c("gene_name", "qglobal_cv_MTLE", "qglobal_cv_MND", 
                                       "wnon_cv_MTLE", "wnon_cv_MND")]
  sig_genes_table <- sig_genes_table[order(pmin(sig_genes_table$qglobal_cv_MTLE, 
                                                sig_genes_table$qglobal_cv_MND)), ]
  print(sig_genes_table)
  write.csv(sig_genes_table, "significant_genes_comparison.csv", row.names = FALSE)
}

if(sum(gene_comparison$sig_either) > 0) {
  # Filter for genes where the significant cohort has non-zero w
  sig_MTLE <- gene_comparison$qglobal_cv_MTLE < 0.05
  sig_MND <- gene_comparison$qglobal_cv_MND < 0.05
  
  sig_with_nonzero_w <- gene_comparison$sig_either & 
    ((sig_MTLE & gene_comparison$wnon_cv_MTLE != 0) |
       (sig_MND & gene_comparison$wnon_cv_MND != 0))
  
  if(sum(sig_with_nonzero_w) > 0) {
    sig_genes_table <- gene_comparison[sig_with_nonzero_w, 
                                       c("gene_name", "qglobal_cv_MTLE", "qglobal_cv_MND", 
                                         "wnon_cv_MTLE", "wnon_cv_MND")]
    
    sig_genes_table <- sig_genes_table[order(pmin(sig_genes_table$qglobal_cv_MTLE, 
                                                  sig_genes_table$qglobal_cv_MND)), ]
    
    print(sig_genes_table)
    write.csv(sig_genes_table, "significant_genes_nonzero_w.csv", row.names = FALSE)
    
    cat(sprintf("\n%d genes with significant q-values and non-zero w in significant cohort\n", 
                sum(sig_with_nonzero_w)))
  } else {
    cat("No significant genes with non-zero w values found\n")
  }
}

dnds_mtle

test <- geneci(dnds_mtle)
test[test$mis_mle > 1 & test$mis_low > 1,]

test[test$tru_mle > 1 & test$tru_low > 1,]

####################
# Plotting
####################

# Plot 1: Global dN/dS comparison
global_long <- melt(global_comparison, id.vars = "group", 
                    variable.name = "mutation_type", value.name = "dNdS")

p1 <- ggplot(global_long, aes(x = mutation_type, y = dNdS, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Global dN/dS Comparison: MTLE vs MND",
       x = "Mutation Type",
       y = "dN/dS ratio",
       fill = "Cohort") +
  theme_bw() +
  theme(text = element_text(size = 12))

p1

ggsave("plot1_global_dnds_comparison.png", p1, width = 8, height = 6, dpi = 300)

# Plot 2: Scatter plot of wnon (nonsense/essential splice) comparison
p2 <- ggplot(gene_comparison, aes(x = wnon_cv_MTLE, y = wnon_cv_MND)) +
  geom_point(aes(color = sig_either), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "blue"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "Significant (q<0.1)")) +
  labs(title = "Gene-level dN/dS Comparison (Nonsense)",
       x = "ωnon (MTLE)",
       y = "ωnon (MND)",
       color = "") +
  theme_bw() +
  theme(text = element_text(size = 12))

p2

ggsave("plot2_wnon_scatter.png", p2, width = 8, height = 7, dpi = 300)

# Plot 3: Scatter plot of wmis (missense) comparison
p3 <- ggplot(gene_comparison, aes(x = wmis_cv_MTLE, y = wmis_cv_MND)) +
  geom_point(aes(color = sig_either), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "blue"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "Significant (q<0.1)")) +
  labs(title = "Gene-level dN/dS Comparison (Missense)",
       x = "ωmis (MTLE)",
       y = "ωmis (MND)",
       color = "") +
  theme_bw() +
  theme(text = element_text(size = 12))

p3

ggsave("plot3_wmis_scatter.png", p3, width = 8, height = 7, dpi = 300)

# Plot 4: Volcano-style plot (delta wnon vs significance)
gene_comparison$min_qglobal <- pmin(gene_comparison$qglobal_cv_MTLE, 
                                    gene_comparison$qglobal_cv_MND)
gene_comparison$neg_log10_q <- -log10(gene_comparison$min_qglobal + 1e-100)

p4 <- ggplot(gene_comparison, aes(x = delta_wnon, y = neg_log10_q)) +
  geom_point(aes(color = sig_either), alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "blue"),
                     labels = c("FALSE" = "Not significant", "TRUE" = "Significant")) +
  labs(title = "Differential Selection: MTLE vs MND",
       x = "Δωnon (MTLE - MND)",
       y = "-log10(min q-value)",
       color = "") +
  theme_bw() +
  theme(text = element_text(size = 12))

p4

ggsave("plot4_differential_selection.png", p4, width = 9, height = 7, dpi = 300)

# Plot 5: Top genes in each cohort (bar plot)
top_n <- 15
top_mtle <- sel_mtle[order(sel_mtle$qglobal_cv), ][1:min(top_n, nrow(sel_mtle)), ]
top_mnd <- sel_mnd[order(sel_mnd$qglobal_cv), ][1:min(top_n, nrow(sel_mnd)), ]

top_combined <- rbind(
  data.frame(gene_name = top_mtle$gene_name, qvalue = top_mtle$qglobal_cv, group = "MTLE"),
  data.frame(gene_name = top_mnd$gene_name, qvalue = top_mnd$qglobal_cv, group = "MND")
)

p5 <- ggplot(top_combined, aes(x = reorder(gene_name, -qvalue), y = -log10(qvalue + 1e-100), fill = group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~group, scales = "free_x", ncol = 1) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") +
  labs(title = paste0("Top ", top_n, " Genes by Selection Significance"),
       x = "Gene",
       y = "-log10(q-value)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 10),
        legend.position = "none")
p5

ggsave("plot5_top_genes.png", p5, width = 12, height = 10, dpi = 300)

cat("\n=== Analysis complete ===\n")
cat("Output files generated:\n")
cat("  - global_dnds_comparison.csv\n")
cat("  - gene_level_comparison.csv\n")
cat("  - significant_genes_comparison.csv\n")
cat("  - plot1_global_dnds_comparison.png\n")
cat("  - plot2_wnon_scatter.png\n")
cat("  - plot3_wmis_scatter.png\n")
cat("  - plot4_differential_selection.png\n")
cat("  - plot5_top_genes.png\n")

p1
p2
p3
p4
p5

####

# Extract total depth from all variants in both cohorts
extract_depth_data <- function(sampleslist1, group_name1, sampleslist2, group_name2) {
  depth_data <- data.frame()
  
  # Process first cohort
  for(i in 1:length(sampleslist1)){
    vcf_object <- read.vcfR(sampleslist1[i], verbose = FALSE)
    
    # Extract sample ID from filename
    SampleID <- basename(sampleslist1[i])
    print(paste0("Processing ", group_name1, " sample: ", SampleID))
    
    # Extract AD field
    gt <- extract.gt(vcf_object, element = "AD")
    ad_split <- strsplit(gt[,1], ",")
    total_depth <- sapply(ad_split, function(x) sum(as.numeric(x)))
    
    # Remove NA values
    total_depth <- total_depth[!is.na(total_depth)]
    
    if(length(total_depth) > 0) {
      depth_sample <- data.frame(
        sampleID = SampleID,
        group = group_name1,
        total_depth = total_depth
      )
      depth_data <- rbind(depth_data, depth_sample)
    }
  }
  
  # Process second cohort
  for(i in 1:length(sampleslist2)){
    vcf_object <- read.vcfR(sampleslist2[i], verbose = FALSE)
    
    # Extract sample ID from filename
    SampleID <- basename(sampleslist2[i])
    print(paste0("Processing ", group_name2, " sample: ", SampleID))
    
    # Extract AD field
    gt <- extract.gt(vcf_object, element = "AD")
    ad_split <- strsplit(gt[,1], ",")
    total_depth <- sapply(ad_split, function(x) sum(as.numeric(x)))
    
    # Remove NA values
    total_depth <- total_depth[!is.na(total_depth)]
    
    if(length(total_depth) > 0) {
      depth_sample <- data.frame(
        sampleID = SampleID,
        group = group_name2,
        total_depth = total_depth
      )
      depth_data <- rbind(depth_data, depth_sample)
    }
  }
  
  return(depth_data)
}

head(depth_data)

depth_summary_table <- depth_data

head(depth_data)

summary_stats <- depth_data %>%
  group_by(sampleID, group) %>%
  summarise(
    n_variants = n(),
    median_depth = median(total_depth),
    mean_depth = mean(total_depth),
    sd_depth = sd(total_depth)
  )

head(summary_stats)

library(dplyr)
library(readr)
library(stringr)

# Define the directory path
dir_path <- "/vast/scratch/users/reid.j/austin_panel/December_Download/nf_run/mosdepth_results_v2/"

# Get the list of bed files
bed_files <- list.files(dir_path, pattern = "\\.per-base\\.bed$", full.names = TRUE)

# Extract sample prefix from filename: only the initial numeric part before dash 
extract_sample_prefix <- function(filename) {
  full_prefix <- str_extract(basename(filename), "^[^\\.]+")
  sample_prefix <- str_extract(full_prefix, "^\\d+")
  return(sample_prefix)
}

# Read all files and calculate median depths
depth_data <- tibble()
for (bed_file in bed_files) {
  sample_prefix <- extract_sample_prefix(bed_file)
  depths <- read_tsv(bed_file, col_names = FALSE, show_col_types = FALSE) |> pull(X3)
  bam_median_depth <- median(depths, na.rm = TRUE)
  depth_data <- depth_data |> bind_rows(tibble(sample_prefix = sample_prefix, bam_median_depth = bam_median_depth))
}

# Prepare summary_stats sample prefix with numeric prefix only
summary_stats <- summary_stats |>
  mutate(sample_prefix = str_extract(sampleID, "^\\d+"))

# Join depth data by matching numeric prefixes only
summary_stats <- left_join(summary_stats, depth_data, by = "sample_prefix") |> select(-sample_prefix)

# Check for successful join
print(summary_stats |> select(sampleID, bam_median_depth) |> head())

head(summary_stats)

length(summary_stats$bam_median_depth.y)

summary_stats

depth_summary_table

# Create depth comparison plot matching age_plot style
plot_depth_comparison <- function(depth_summary_table, title = "Total Depth Distribution Comparison") {
  
  # Create boxplot with matching style to age_plot
  p <- ggplot(depth_summary_table, aes(x = group, y = bam_median_depth.y, color = group)) +
    geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
    scale_y_continuous(limits = c(0, 1000), breaks = seq(0, 1000, by = 100)) +
    labs(title = title,
         x = "Cohort",
         y = "Median depth",
         color = "Cohort") +
    theme_bw(base_size = 16) +  # Match base_size from age_plot
    theme(
      legend.position = "right",  # Match age_plot
      plot.title = element_text(size = 18, face = "bold"),  # Match age_plot
      axis.title = element_text(size = 16),  # Match age_plot
      axis.text = element_text(size = 14),  # Match age_plot
      legend.title = element_text(size = 16),  # Match age_plot
      legend.text = element_text(size = 14)  # Match age_plot
    )
  
  return(list(plot = p, summary = summary_stats))
}

# Usage example:
result <- plot_depth_comparison(summary_stats, "Median Depth Distribution by Cohort")
depth_plot <- result$plot
print(depth_plot)

colnames(summary_stats)
# Usage example:
# Extract depth data
depth_data <- extract_depth_data(mtle, "MTLE", mnd, "MND")

head(depth_data)

# Create violin plot and get summary statistics
result <- plot_depth_comparison(depth_data)
result$plot

# Save plot
ggsave("depth_comparison_violin.png", result$plot, width = 8, height = 6, dpi = 300)

# Perform statistical test
wilcox_test <- wilcox.test(total_depth ~ group, data = depth_data)
print("Wilcoxon rank-sum test:")
print(wilcox_test)

