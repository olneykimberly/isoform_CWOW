# Load the necessary libraries
library(tidyverse)

# 1. Read the data
isoform_data <- read_tsv("/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/collapsed_isos/NA18-285_LBD_collapsed_sorted_classification.filtered_lite_classification.txt")
test <- read_tsv("/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/collapsed_isos/NA15-031_CONTROL_collapsed_sorted_classification.filtered_lite_classification.txt")

# Define the list of sample IDs
sampleIDs <- c('NA15-333_CONTROL', 'NA20-121_CONTROL', 'NA15-031_CONTROL', 'NA18-285_LBD', 'NA19-290_LBD', 'NA07-150_LBD') 

# Define the base path to your files
base_path <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/collapsed_isos/"

# Create an empty data frame to store the results
counts_df <- data.frame(sample_id = character(), 
                        line_count = numeric(), 
                        stringsAsFactors = FALSE)

# Use a for loop to iterate through each sample ID
for (sample in sampleIDs) {
  # Construct the full file path
  file_path <- paste0(base_path, sample, "_collapsed_sorted_classification.filtered_lite_classification.txt")
  
  # Check if the file exists before attempting to read it
  if (file.exists(file_path)) {
    # Read the file and count the number of lines
    # Using read_lines() is more efficient than read_tsv() for just counting lines
    lines <- read_lines(file_path)
    line_count <- length(lines)
    
    # Add the result to our data frame
    counts_df <- bind_rows(counts_df, 
                           data.frame(sample_id = sample, line_count = line_count))
  } else {
    # Print a message if the file is not found
    message(paste("Warning: File not found for sample:", sample))
  }
}

# 4. Create the bar plot
# Ensure counts_df is not empty before plotting
if (nrow(counts_df) > 0) {
  ggplot(counts_df, aes(x = sample_id, y = line_count)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Number of Isoforms Per Sample",
      x = "Sample",
      y = "Number of Isoforms (Lines)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
} else {
  message("No data to plot. Please check your sample IDs and file paths.")
}

table(isoform_data$structural_category)
table(isoform_data$subcategory)
table(isoform_data$coding)