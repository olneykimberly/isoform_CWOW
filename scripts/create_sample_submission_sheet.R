library(dplyr)
# Sample submission information
samples_selected <- read.delim("metadata/long_read_sample_selection.txt")

# Freezer information
freezer_key <- read.delim("metadata/freezer_key.tsv")
sample_freezer <- read.delim("metadata/CWOW_plate_well_RIN.tsv")

# merge
samples_selected_freezer_info <- 
  merge(samples_selected, sample_freezer, by = "NPID")

metadata <- merge(
  x = samples_selected_freezer_info,
  y = freezer_key[, c("Freezer_rack", "Freezer_box", "CWOW_batch")],
  by.x = "Box.x",
  by.y = "CWOW_batch",
  all.x = TRUE
)

rm(samples_selected, freezer_key, sample_freezer, samples_selected_freezer_info)
metadata$NPID_no_dash <- gsub("-", "", metadata$NPID)
metadata$ATS_with_dash <- gsub("_", "-", metadata$group.redefined)
metadata$Study_Subject_ID <- paste0(metadata$NPID_no_dash, "_", metadata$ATS_with_dash)
metadata$Species <- "Homo Sapiens"
metadata$Collection_Date <- "1/1/2001"
metadata$Archive_Type <- "Isolation"
metadata$Study_Specimen_ID <- metadata$NPID_no_dash
metadata$Specimen_Source <- "Brain"
metadata$Specimen_Subgroup <- "Constitutional"
metadata$Specimen_Type <- "RNA"
metadata$Specimen_Fraction <- "Whole"
metadata$Isolation_Method <- "RNA - Qiagen RNeasy Universal plus"
metadata$Isolation_Buffer <- "Water"
metadata$TGen_ID <- ""
metadata$Volume_Units <- "ul"
metadata$Tapestation_DIN_or_RIN <- metadata$RIN.x
# sex -> Gender
metadata$Gender <- metadata$sex_inferred
metadata$Tapestation_Conc <- metadata$Tapestation_Conc
metadata$Internal_Tapestation_DV200
#  Race -> Ethnicity
metadata$Ethnicity <- metadata$Race 
metadata$Refrigerator <- ""
metadata$Rack <- ""
metadata$Box <- ""
metadata$Slot <- ""
metadata$Nanodrop_Conc <- ""
metadata$Nanodrop_260_280 <- ""
metadata$Nanodrop_260_230 <- ""
metadata$Qubit_Assay <- ""
metadata$Qubit_Conc <- ""
metadata$Tapestation_Assay <- ""

metadata$Notes <- ""
metadata$Import_Level_CSC_USE_ONLY <- ""
metadata$batch_box <- metadata$Box.y
metadata$group <- metadata$group.redefined_again
sample_submission_info <- metadata[, c("Study_Subject_ID", "Species", "Gender", "Ethnicity", "Collection_Date", 
                                       "Archive_Type", "Study_Specimen_ID", "TGen_ID", "Specimen_Source", 
                                       "Specimen_Subgroup", "Specimen_Type", "Specimen_Fraction", "Isolation_Method",
                                       "Isolation_Buffer", "Volume", "Volume_Units", 
                                       "Refrigerator", "Rack", "Box", "Slot", "Nanodrop_Conc", "Nanodrop_260_280", "Nanodrop_260_230", "Qubit_Assay", "Qubit_Conc", "Tapestation_Assay",
                                       "Tapestation_DIN_or_RIN", "Tapestation_Conc", "Internal_Tapestation_DV200", 
                                       "Notes", "Import_Level_CSC_USE_ONLY", 
                                       "Freezer_rack", "Freezer_box", "batch_box", "well_ID", "NPID", "group", "Braak.NFT", "Thal.amyloid", "Cing.LB", "VaD", "TDP.43", "APOE", "Brain.wt", "Age", "PMI")]


# Already sequenced to 20Million reads
# 'NA07-150_LBD', 'NA15-031_CONTROL', 'NA15-333_CONTROL', 'NA19-290_LBD', 'NA18-285_LBD', 'NA20-121_CONTROL'

set.seed(28) 

samples_to_exclude <- sample_submission_info %>%
  filter(group == "LBD_ATS", Gender == "male") %>%
  sample_n(9) %>%
  pull("NPID")

final_data <- sample_submission_info %>%
  filter(!NPID %in% samples_to_exclude)

write.table(final_data, "metadata/sample_submission_info.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Create a key for group and sex 
final_data <- final_data %>%
  mutate(Group_Sex = paste(group, Gender, sep = "_"))

df_batched <- final_data %>%
  mutate(
    row_num = row_number(), # Global sequential numbering after shuffle
    SMRTcell_ID = as.character(ceiling(row_num / 6))
  ) %>%
  # --- Add SMRTcell_sample_number ---
  # Group by the newly assigned SMRTcell_ID
  group_by(SMRTcell_ID) %>%
  # Assign a sequential number (1 to 6) within each SMRTcell_ID group
  mutate(
    SMRTcell_sample_number = row_number()
  ) %>%
  ungroup() 

write.table(df_batched, "metadata/sample_submission_info_with_SMRTcell.txt", sep = "\t", quote = FALSE, row.names = FALSE)

SMRT_sex <- as.data.frame(table(df_batched$SMRTcell_ID, df_batched$Gender))
SMRT_group <- as.data.frame(table(df_batched$SMRTcell_ID, df_batched$group))

#----------
library(dplyr)

samples <- read.delim("/Users/kolney/Desktop/sample_updated_batches.txt")
samples$Sample_number <- NULL
samples$SMRTcell_ID <- NULL
samples$SMRTcell_sample_number <- NULL


df_batched <- samples %>%
  mutate(
    row_num = row_number(), # Global sequential numbering after shuffle
    SMRTcell_ID = as.character(ceiling(row_num / 6))
  ) %>%
  # --- Add SMRTcell_sample_number ---
  # Group by the newly assigned SMRTcell_ID
  group_by(SMRTcell_ID) %>%
  # Assign a sequential number (1 to 6) within each SMRTcell_ID group
  mutate(
    SMRTcell_sample_number = row_number()
  ) %>%
  ungroup() 

SMRT_sex <- as.data.frame(table(df_batched$SMRTcell_ID, df_batched$Gender))
SMRT_group <- as.data.frame(table(df_batched$SMRTcell_ID, df_batched$group))
ID <- as.data.frame(table(df_batched$Study_Subject_ID))


write.table(df_batched, "/Users/kolney/Desktop/sample_updated_batches_updated.txt", sep = "\t", quote = FALSE, row.names = FALSE)

