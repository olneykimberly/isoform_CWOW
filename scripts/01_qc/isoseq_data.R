# Load the necessary libraries
# If you don't have them installed, run: install.packages("tidyverse")
library(tidyverse)

# Define the list of sample IDs
# Replace these with your actual sample IDs
sampleIDs <- c('NA15-333_CONTROL', 'NA20-121_CONTROL', 'NA15-031_CONTROL', 'NA18-285_LBD', 'NA19-290_LBD', 'NA07-150_LBD') 

# Define the base path to your files
base_path <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/collapsed_isos/"

# Define the columns of interest for which to create frequency tables
category_cols <- c("structural_category", "subcategory", "coding")

# Create empty lists to store the tables for each category
category_tables <- list()
subcategory_tables <- list()
coding_tables <- list()

# Create a master data frame to store classification counts for plotting
# The columns will be: sample_id, category, feature, count
# Initialize as an empty tibble for proper Tidyverse compatibility
classification_counts_df <- tibble(sample_id = character(), category = character(), feature = character(), count = numeric())

# Create a data frame for non_coding counts specifically
non_coding_df <- tibble(sample_id = character(), non_coding_count = numeric())

# Use a for loop to iterate through each sample ID
for (sample in sampleIDs) {
  # Construct the full file path
  file_path <- paste0(base_path, sample, "_collapsed_sorted_classification.filtered_lite_classification.txt")
  
  # Check if the file exists before attempting to read it
  if (file.exists(file_path)) {
    
    # Read the isoform classification data (assuming it's a tab-separated file)
    isoform_data <- tryCatch(
      read_tsv(file_path, show_col_types = FALSE),
      error = function(e) {
        message(paste("Error reading file for sample:", sample, "- ", e$message))
        return(NULL)
      }
    )
    
    if (!is.null(isoform_data)) {
      
      # 1. Create and store a table of structural_category
      struct_table <- as.data.frame(table(isoform_data$structural_category))
      colnames(struct_table) <- c("structural_category", "Freq")
      category_tables[[sample]] <- struct_table
      message(paste("Structural Category Table for sample:", sample))
      print(struct_table)
      
      # Prepare structural_category data for plotting (Appended to master DF)
      temp_struct_df <- struct_table %>%
        transmute(
          sample_id = sample,
          category = "structural_category",
          feature = structural_category,
          count = Freq
        )
      classification_counts_df <- bind_rows(classification_counts_df, temp_struct_df)
      
      # 2. Create and store a table of subcategory
      subcat_table <- as.data.frame(table(isoform_data$subcategory))
      colnames(subcat_table) <- c("subcategory", "Freq")
      subcategory_tables[[sample]] <- subcat_table
      message(paste("Subcategory Table for sample:", sample))
      print(subcat_table)
      
      # **Prepare subcategory data for plotting (Appended to master DF)**
      temp_subcat_df <- subcat_table %>%
        transmute(
          sample_id = sample,
          category = "subcategory",
          feature = subcategory,
          count = Freq
        )
      classification_counts_df <- bind_rows(classification_counts_df, temp_subcat_df)
      
      # 3. Create and store a table of coding
      coding_table <- as.data.frame(table(isoform_data$coding))
      colnames(coding_table) <- c("coding", "Freq")
      coding_tables[[sample]] <- coding_table
      message(paste("Coding Table for sample:", sample))
      print(coding_table)
      
      # **Extract non_coding count for separate plot**
      non_coding_count <- coding_table %>%
        filter(coding == "non_coding") %>%
        pull(Freq)
      
      # If 'non_coding' exists, store the count, otherwise store 0
      non_coding_count <- ifelse(length(non_coding_count) > 0, non_coding_count, 0)
      
      non_coding_df <- bind_rows(non_coding_df, 
                                 tibble(sample_id = sample, non_coding_count = non_coding_count))
      
    }
    
  } else {
    # Print a message if the file is not found
    message(paste("Warning: File not found for sample:", sample))
  }
}

# Define a placeholder for project_ID and ensure the results directory exists
project_ID <- "isoseq"
results_dir <- paste0("results_", project_ID)
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# ----------------------------------------------------------------------
## Relative Abundance Plot for Subcategory
# ----------------------------------------------------------------------

if (nrow(classification_counts_df) > 0) {
  
  # Calculate relative abundance for subcategory
  subcategory_plot_data <- classification_counts_df %>%
    filter(category == "subcategory") %>%
    group_by(sample_id) %>%
    # Calculate proportion within each sample
    mutate(relative_abundance = count / sum(count)) %>%
    ungroup()
  
  # Save the plot
  pdf(paste0(results_dir, "/subcategory_relative_abundance.pdf"), width = 10, height = 7)
  
  subcategory_plot <- ggplot(subcategory_plot_data, aes(x = sample_id, y = relative_abundance, fill = feature)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    scale_y_continuous(labels = scales::percent) + 
    theme_minimal(base_size = 14) +
    labs(
      title = "Relative Abundance of Isoform Subcategories",
      x = "Sample ID",
      y = "Relative Abundance",
      fill = "Subcategory"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  print(subcategory_plot)
  dev.off()
  
  message("Relative abundance plot for subcategory saved to: ", paste0(results_dir, "/subcategory_relative_abundance.pdf"))
  
} else {
  message("No data available to plot subcategory relative abundance.")
}

# ----------------------------------------------------------------------
## Bar Plot of Non-Coding Isoform Counts
# ----------------------------------------------------------------------

if (nrow(non_coding_df) > 0) {
  
  # Sort data for better visualization
  non_coding_df <- non_coding_df %>%
    mutate(sample_id = reorder(sample_id, -non_coding_count))
  
  # Save the plot
  pdf(paste0(results_dir, "/non_coding_counts_per_sample.pdf"), width = 6, height = 6)
  
  non_coding_plot <- ggplot(non_coding_df, aes(x = sample_id, y = non_coding_count)) +
    geom_bar(stat = "identity", fill = "salmon", color = "black") +
    geom_text(aes(label = non_coding_count), vjust = -0.5, size = 4) + # Add count labels
    theme_minimal(base_size = 14) +
    labs(
      title = "Count of Non-Coding Isoforms",
      x = "Sample ID",
      y = "Number of Non-Coding Isoforms"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(non_coding_plot)
  dev.off()
  
  message("Bar plot for non-coding counts saved to: ", paste0(results_dir, "/non_coding_counts_per_sample.pdf"))
  
} else {
  message("No data available to plot non-coding counts.")
}

# ----------------------------------------------------------------------
## (Optional) Structural Category Plot (kept from previous iteration)
# ----------------------------------------------------------------------

# NOTE: The code for the structural_category plot remains from the previous iteration 
# and is not strictly needed in this final output unless requested, but is left in the
# master DF creation for completeness. To run the structural category plot, you'd 
# reuse the plotting logic from the previous prompt, filtering for 'structural_category'.

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
pdf(paste0("results_", project_ID, "/isoforms_per_sample.pdf"), width = 6, height = 6)
if (nrow(counts_df) > 0) {
  ggplot(counts_df, aes(x = sample_id, y = line_count)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Number of Isoforms",
      x = "Sample",
      y = "Number of Isoforms (Lines)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
} else {
  message("No data to plot. Please check your sample IDs and file paths.")
}
dev.off()


