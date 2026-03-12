#!/usr/bin/env Rscript

# Define file paths
input_file  <- "/vast/scratch/users/reid.j/duplex_sequencing_workflow/mbc002_ids.tsv"
output_file <- "/vast/scratch/users/reid.j/duplex_sequencing_workflow/mbc002_ids_plus_vcfpath_plus_bam.tsv"
search_dir  <- "/vast/scratch/users/reid.j/austin_panel"

# Read the input TSV file — no header, two columns
cat("Reading input file:", input_file, "\n")
sample_data <- read.delim(input_file, header = FALSE, stringsAsFactors = FALSE, col.names = c("Sample_ID", "Group"))

cat("Samples found:", nrow(sample_data), "\n")
print(sample_data)

# Initialize the Filepath columns with NA
sample_data$VCF_Filepath <- NA
sample_data$BAM_Filepath <- NA

# Helper: read a line from stdin, handling EOF gracefully
read_stdin <- function() {
  line <- readLines(con = stdin(), n = 1)
  if (length(line) == 0) return("")
  return(trimws(line))
}

# Helper: validated numeric choice from user
get_user_choice <- function(prompt_msg, max_choice) {
  repeat {
    cat(prompt_msg)
    input      <- read_stdin()
    choice_num <- suppressWarnings(as.integer(input))
    if (length(choice_num) == 1 && !is.na(choice_num) &&
        choice_num >= 1 && choice_num <= max_choice) {
      return(choice_num)
    }
    cat("Invalid choice. Please enter a number between 1 and", max_choice, "\n")
  }
}

# Iterate over each row
for (i in 1:nrow(sample_data)) {
  sample_id <- sample_data$Sample_ID[i]
  group     <- sample_data$Group[i]

  cat("\n----------------------------------------\n")
  cat("Processing row", i, "of", nrow(sample_data), "\n")
  cat("Sample ID:", sample_id, "\n")
  cat("Group:", group, "\n")

  if (is.na(sample_id) || nchar(trimws(sample_id)) == 0) {
    cat("WARNING: Empty sample ID at row", i, "- skipping\n")
    next
  }

  # Search for matching VCF files using full sample ID (e.g. "51146-2")
  vcf_pattern    <- paste0("^", sample_id, ".*\\.hard-filtered\\.vcf\\.gz$")
  matching_files <- list.files(
    path       = search_dir,
    pattern    = vcf_pattern,
    full.names = TRUE,
    recursive  = TRUE
  )

  if (length(matching_files) == 0) {
    cat("WARNING: No matching VCF files found for", sample_id, "\n")
    next
  } else if (length(matching_files) == 1) {
    cat("Found 1 matching VCF file:\n  ", matching_files, "\n")
    sample_data$VCF_Filepath[i] <- matching_files
  } else {
    cat("Found", length(matching_files), "matching VCF files:\n")
    for (j in seq_along(matching_files)) cat("  [", j, "] ", matching_files[j], "\n", sep = "")
    idx <- get_user_choice(
      paste0("Select VCF file for ", sample_id, " (1-", length(matching_files), "): "),
      length(matching_files)
    )
    sample_data$VCF_Filepath[i] <- matching_files[idx]
    cat("Selected:", matching_files[idx], "\n")
  }

  # Search for BAM in the same directory as the selected VCF
  if (!is.na(sample_data$VCF_Filepath[i])) {
    vcf_dir       <- dirname(sample_data$VCF_Filepath[i])
    cat("\nSearching for BAM file in:", vcf_dir, "\n")

    bam_pattern   <- paste0("^", sample_id, ".*\\.bam$")
    matching_bams <- list.files(
      path       = vcf_dir,
      pattern    = bam_pattern,
      full.names = TRUE,
      recursive  = FALSE
    )

    if (length(matching_bams) == 0) {
      cat("WARNING: No matching BAM files found in", vcf_dir, "\n")
    } else if (length(matching_bams) == 1) {
      cat("Found 1 matching BAM file:\n  ", matching_bams, "\n")
      sample_data$BAM_Filepath[i] <- matching_bams
    } else {
      cat("Found", length(matching_bams), "matching BAM files:\n")
      for (j in seq_along(matching_bams)) cat("  [", j, "] ", matching_bams[j], "\n", sep = "")
      idx <- get_user_choice(
        paste0("Select BAM file for ", sample_id, " (1-", length(matching_bams), "): "),
        length(matching_bams)
      )
      sample_data$BAM_Filepath[i] <- matching_bams[idx]
      cat("Selected:", matching_bams[idx], "\n")
    }
  }
}

# Summary
cat("\n========================================\n")
cat("Processing complete!\n")
cat("Total samples processed:", nrow(sample_data), "\n")
cat("VCF files found:        ", sum(!is.na(sample_data$VCF_Filepath)), "\n")
cat("BAM files found:        ", sum(!is.na(sample_data$BAM_Filepath)), "\n")
cat("VCF files not found:    ", sum(is.na(sample_data$VCF_Filepath)),  "\n")
cat("BAM files not found:    ", sum(is.na(sample_data$BAM_Filepath)),  "\n")

# Write output
cat("\nWriting output to:", output_file, "\n")
write.table(sample_data, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Done!\n\n")
print(head(sample_data))