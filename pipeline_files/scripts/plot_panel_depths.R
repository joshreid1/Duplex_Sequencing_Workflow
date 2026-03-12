# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Define the directories containing the files
consensus_dir <- "/vast/scratch/users/reid.j/austin_panel/nextflow_run/mosdepth_results"
raw_dir <- ""

# List all files in the directories matching the pattern
file_pattern <- "\\.gene-list\\.per-base\\.bed$"
consensus_files <- list.files(consensus_dir, pattern = file_pattern, full.names = TRUE)
raw_files <- list.files(raw_dir, pattern = file_pattern, full.names = TRUE)

# Helper function to read and format a file, including the Processing variable
read_and_format <- function(file_path, processing_type) {
  # Extract sample ID, repeat type, and version
  file_name <- basename(file_path)
  sample_id <- strsplit(file_name, "\\.")[[1]][1]
  batch_id <- strsplit(file_name, "\\.")[[1]][2]
  repeat_type <- sub(".*\\.(FRESH|FFPE)\\..*", "\\1", file_name) # Extract repeat type (FRESH or FFPE)
  
  # Read the data
  df <- read.table(file_path, sep = "\t", header = FALSE, 
                   col.names = c("CHROM", "POS", "COVERAGE"))
  
  # Add sample ID, repeat type, and processing type columns
  df <- df %>%
    mutate(Sample_ID = sample_id, Batch_ID = batch_id, Repeat = repeat_type, Processing = processing_type)
  
  return(df)
}

# Read and combine data from consensus directory
consensus_df <- bind_rows(lapply(consensus_files, read_and_format, "Consensus"))

summary_df <- consensus_df %>%
  group_by(Sample_ID, Batch_ID, Repeat, Processing) %>%
  summarise(
    coverage_mean = mean(COVERAGE, na.rm = TRUE),
    coverage_median = median(COVERAGE, na.rm = TRUE),
    coverage_min = min(COVERAGE, na.rm = TRUE),
    coverage_max = max(COVERAGE, na.rm = TRUE),
    coverage_Q1 = quantile(COVERAGE, 0.25, na.rm = TRUE),
    coverage_Q3 = quantile(COVERAGE, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )

# View result
head(summary_df)

sort(summary_df$coverage_median)

median(summary_df$coverage_median)

write.table(summary_df, file = "coverage_summary.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)


# Check if raw files exist and if not, only use consensus data
if (length(raw_files) > 0) {
  raw_df <- bind_rows(lapply(raw_files, read_and_format, "Raw"))
  combined_df <- bind_rows(consensus_df, raw_df)
} else {
  combined_df <- consensus_df
}

combined_df$Repeat <- 'Fresh'

# Ensure the Repeat column has the correct order
combined_df <- combined_df %>%
  mutate(Repeat = factor(Repeat, levels = c("Fresh", "FFPE", "FFPE-AustinPath")))

#combined_df[combined_df$Sample_ID == "HORIZON", "Repeat"] = "FFPE-AustinPath"
#combined_df[combined_df$Sample_ID == "M23017", "Repeat"] = "FFPE-AustinPath"

plot_df <- combined_df

# Filter for the specific sample IDs
#plot_df <- combined_df %>%
#  filter(Sample_ID %in% c('37673-6', '37977-4', '38229-4', 'HORIZON'))


# Ensure Processing is ordered as Raw first (if raw data exists)
if (length(raw_files) > 0) {
  plot_df$Processing <- factor(plot_df$Processing, levels = c("Raw", "Consensus"))
} else {
  plot_df$Processing <- factor(plot_df$Processing, levels = c("Consensus", "Raw"))
}


#plot_df <- plot_df %>%
#  mutate(Sample_ID = factor(Sample_ID, levels = c("Sample 1", "Sample 2", "Sample 3", "Horizon Control")))

# Calculate median coverage for each Sample_ID and Processing combination
median_order <- plot_df %>%
  group_by(Sample_ID, Processing) %>%
  summarise(median_coverage = median(COVERAGE, na.rm = TRUE)) %>%
  arrange(median_coverage) %>%
  pull(Sample_ID) %>%
  unique()

# Order Sample_ID factor based on median coverage
plot_df$Sample_ID <- factor(plot_df$Sample_ID, levels = median_order)

# Create a combined label for strip text
batch_labels <- plot_df %>%
  dplyr::distinct(Batch_ID) %>%
  dplyr::mutate(
    Batch_Label = paste0(
      Batch_ID, "\n"
    )
  )

# Merge back to main data for faceting
plot_df <- plot_df %>%
  dplyr::left_join(batch_labels, by = "Batch_ID")

boxplot <- plot_df %>%
  filter(!is.na(Batch_Label) & Batch_Label != "") %>%
  ggplot(aes(x = factor(Sample_ID), y = COVERAGE, fill = Processing)) +
  geom_boxplot(width = 0.4, position = position_dodge(width = 0.6), outlier.shape = NA) +
  theme_minimal() +
  scale_y_continuous(
    limits = c(0, 2000),
    breaks = scales::pretty_breaks(n = 10),
    sec.axis = dup_axis(name = NULL)
  ) +
  scale_fill_manual(
    values = c("Consensus" = "#e4980f"),
    labels = c("Consensus reads")
  ) +
  labs(y = "Depth of Coverage") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    panel.spacing.x = unit(1, "lines"),
    panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    # Remove the box and fully display labels
    strip.text.x = element_text(
      size = 10,
      face = "bold",
      margin = margin(t = 5, b = 5, l = 0, r = 0, unit = "pt")
    ),
    strip.background = element_blank(), # This removes the background box
    # Adjust plot margins
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
  ) +
  facet_grid(. ~ Batch_Label, scales = "free_x", space = "free_x") +
  coord_cartesian(clip = "off")

# Enhanced width calculation
batch_labels <- unique(plot_df$Batch_Label[!is.na(plot_df$Batch_Label) & plot_df$Batch_Label != ""])
n_batches <- length(batch_labels)
max_label_length <- max(nchar(batch_labels), na.rm = TRUE)

# More generous width calculation
min_width_per_batch <- pmax(6, max_label_length * 0.2)  # Increased base width and label factor
total_width <- max(24, n_batches * min_width_per_batch)  # Higher minimum

print(boxplot)

# Save as PNG
ggsave(
  filename = "coverage_boxplot.png",
  plot = boxplot,
  width = total_width,     # adjust width (in inches)
  height = 8,     # adjust height (in inches)
  units = "in",
  dpi = 300       # publication-quality resolution
)

# Create a violin plot
violin_plot <- ggplot(plot_df, aes(x = Sample_ID, y = COVERAGE, fill = Processing)) +
  geom_violin(width = 0.8, position = position_dodge(width = 0.8), trim = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.6) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 4000), breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values = c("Raw" = "#4472c4", "Consensus" = "#e4980f"),
                    labels = c("Raw mapped reads", "Consensus reads")) +
  labs(y = "Depth of Coverage") + # No x-axis label or legend title
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    panel.spacing.x = unit(1, "lines"),
    panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),
    plot.background = element_rect(color = "#4472c4", fill = NA, linewidth = 1.5),
    legend.position = "right", # Move the legend to the right
    legend.title = element_blank(), # Remove legend title
    legend.text = element_text(size = 14) # Increase font size for legend labels
  )

# Print the violin plot
print(violin_plot)