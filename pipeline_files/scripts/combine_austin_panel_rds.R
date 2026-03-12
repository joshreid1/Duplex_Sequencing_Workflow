#!/usr/bin/env Rscript
# combine_rds_to_excel.R
# This script iterates over all RDS files in a specified directory,
# creates an Excel workbook with one worksheet per sample plus a summary sheet,
# and writes the output workbook to an .xlsx file.

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

# Check and collect command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript combine_rds_to_excel.R <input_directory> <output_excel_file>")
}

input_dir <- args[1]
output_excel <- args[2]

# List all RDS files in the directory
rds_files <- list.files(path = input_dir, pattern = "\\.RDS$", full.names = TRUE)
if(length(rds_files) == 0){
  stop("No RDS files found in the specified directory.")
}

# Create an Excel workbook
wb <- createWorkbook()

# Function to add a worksheet with a unique sheet name
addUniqueWorksheet <- function(wb, sheet_name, data) {
  base_name <- sheet_name
  suffix <- 1
  # Ensure the worksheet name is unique (openxlsx will complain if duplicate)
  while(base_name %in% names(wb)){
    base_name <- paste0(sheet_name, "_", suffix)
    suffix <- suffix + 1
  }

  # Truncate to 31 characters if too long
  if(nchar(base_name) > 31) {
    base_name <- substr(base_name, 1, 31)
  }

  addWorksheet(wb, base_name)
  writeDataTable(wb, sheet = base_name, x = data, tableStyle = "TableStyleLight9")
}

# Initialize list for summary data
all_variants <- list()

# Process each RDS file
for(rds_file in rds_files){
  message("Processing file: ", rds_file)
  # Read the RDS file (assumed to be a data.frame)
  sample_data <- readRDS(rds_file)
  
  # Extract a sample name from the file name (without extension)
  sample_name <- tools::file_path_sans_ext(basename(rds_file))
  
  # Handle empty data frames
  if(nrow(sample_data) == 0){
    message("  Warning: ", sample_name, " is empty. Creating placeholder sheet.")
    # Create a minimal data frame with one row indicating no data
    sample_data <- data.frame(
      Sample_ID = sample_name,
      Note = "No variants found",
      stringsAsFactors = FALSE
    )
  } else {
    # Add Sample_ID column only if data exists and column is missing
    if(!"Sample_ID" %in% colnames(sample_data)){
      sample_data$Sample_ID <- sample_name
    }
  }
  
  # Add the sample's data as a new worksheet in the workbook
  addUniqueWorksheet(wb, sample_name, sample_data)
  
  # Append the data to the summary list (only if not a placeholder)
  if(!"Note" %in% colnames(sample_data)){
    all_variants[[sample_name]] <- sample_data
  }
}

# Combine all data frames into one summary data frame
combined_variants <- bind_rows(all_variants)

# Add the summary sheet to the workbook
addWorksheet(wb, "Summary")
writeDataTable(wb, sheet = "Summary", x = combined_variants, tableStyle = "TableStyleLight9")

# Save the workbook to the specified Excel file
saveWorkbook(wb, output_excel, overwrite = TRUE)
cat("Excel workbook created successfully at", output_excel, "\n")
