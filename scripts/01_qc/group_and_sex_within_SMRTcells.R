# Load the necessary library
library(ggplot2)

# --- 1. Load Data ---
# Define the file path
file_path <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/metadata/finale_n288.tsv"

# Read the data into a data frame
# Assuming the file is tab-separated as implied by the .tsv extension
# and has a header row.
metadata <- read.delim(file_path, header = TRUE, sep = "\t")
table(metadata$Gender)
# --- 2. Data Preparation/Cleanup (Optional but Recommended) ---
# Ensure the grouping columns are factors for proper plotting
# The 'SMRTcell_GROUPS' column for the x-axis (24 levels)
# The 'group' column for the fill color (5 levels)
# The 'Gender' column for the second plot's fill color (2 levels)
metadata$SMRTcell_GROUPS <- as.factor(metadata$SMRTcell_GROUPS)
metadata$group <- as.factor(metadata$group)
metadata$Gender <- as.factor(metadata$Gender)

# Ensure RIN_new is treated as a numeric variable
# RIN stands for RNA Integrity Number, which is a continuous measure.
metadata$RIN_new <- as.numeric(metadata$RIN_new)
metadata$Age <- as.numeric(metadata$Age)


# =======================================================
## 🎻 Violin Plot of RIN_new across SMRTcell_GROUPS
# =======================================================

# =======================================================
## 📈 Box Plot and Jitter Plot of RIN_new across SMRTcell_GROUPS
# =======================================================

plot_rin_box_jitter <- ggplot(metadata, aes(x = SMRTcell_GROUPS, y = RIN_new)) +
  
  # 1. Add the Box Plot (summary statistics)
  # 'outlier.shape = NA' prevents the box plot from plotting outliers, 
  # as the jitter plot will show all points.
  geom_boxplot(
    fill = "lightgray", 
    color = "black", 
    alpha = 0.7, 
    outlier.shape = NA
  ) +
  
  # 2. Add the Individual Data Points (raw data)
  # 'width' controls the spread of points horizontally.
  # 'height' can be left as 0 to avoid vertical jitter, which would distort RIN values.
  geom_jitter(
    aes(color = SMRTcell_GROUPS), # Color by SMRTcell Group for differentiation
    width = 0.2,                  # Controls horizontal spread
    size = 1.5,                   # Controls point size
    alpha = 0.8
  ) +
  
  # Add labels and title
  labs(
    title = "RIN_new Distribution by SMRTcell Group",
    x = "SMRTcell Group",
    y = "RIN_new Value (0-10)",
    color = "SMRTcell Group"
  ) +
  
  # Enforce the known RIN range on the y-axis (0 to 10)
  ylim(0, 10) +
  
  # Improve aesthetics: Rotate x-axis labels, remove the redundant legend
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none" # Hide the color legend
  ) +
  
  # Use a suitable color palette for the 24 groups
  scale_color_viridis_d() # Requires 'viridis' package to be loaded, or use 'Paired'

# You might need to install the 'viridis' package if you use this scale: install.packages("viridis")
# If you don't want to install 'viridis', you can use the built-in:
# scale_color_brewer(palette = "Paired") 

# Display the box and jitter plot
print(plot_rin_box_jitter)
#

plot_rin_box_jitter <- ggplot(metadata, aes(x = SMRTcell_GROUPS, y = Age)) +
  
  # 1. Add the Box Plot (summary statistics)
  # 'outlier.shape = NA' prevents the box plot from plotting outliers, 
  # as the jitter plot will show all points.
  geom_boxplot(
    fill = "lightgray", 
    color = "black", 
    alpha = 0.7, 
    outlier.shape = NA
  ) +
  
  # 2. Add the Individual Data Points (raw data)
  # 'width' controls the spread of points horizontally.
  # 'height' can be left as 0 to avoid vertical jitter, which would distort RIN values.
  geom_jitter(
    aes(color = SMRTcell_GROUPS), # Color by SMRTcell Group for differentiation
    width = 0.2,                  # Controls horizontal spread
    size = 1.5,                   # Controls point size
    alpha = 0.8
  ) +
  
  # Add labels and title
  labs(
    title = "Age Distribution by SMRTcell Group",
    x = "SMRTcell Group",
    y = "Age Value (0-10)",
    color = "SMRTcell Group"
  ) +
  
  # Enforce the known RIN range on the y-axis (0 to 10)
  ylim(50, 105) +
  
  # Improve aesthetics: Rotate x-axis labels, remove the redundant legend
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none" # Hide the color legend
  ) +
  
  # Use a suitable color palette for the 24 groups
  scale_color_viridis_d() # Requires 'viridis' package to be loaded, or use 'Paired'

# You might need to install the 'viridis' package if you use this scale: install.packages("viridis")
# If you don't want to install 'viridis', you can use the built-in:
# scale_color_brewer(palette = "Paired") 

# Display the box and jitter plot
print(plot_rin_box_jitter)


# =======================================================
# --- PLOT 1: Stacked Barplot of 'group' within 'SMRTcell_GROUPS' ---
# =======================================================

# 24 bars on the x-axis, 5 colors per bar
plot_group_smrtcell <- ggplot(metadata, aes(x = SMRTcell_GROUPS, fill = group)) +
  geom_bar(position = "stack") +
  
  # Add labels and title
  labs(
    title = "Distribution of 'group' within SMRTcell Groups",
    x = "SMRTcell Group",
    y = "Count of Samples",
    fill = "Group" # Legend title for the fill color
  ) +
  
  # Improve aesthetics: rotate x-axis labels for readability 
  # given 24 levels and a simple theme
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  
  # Optional: Use a specific color palette (e.g., 'Set1' which has 9 colors)
  # or 'Paired' if you want a different look
  scale_fill_brewer(palette = "Set1")

# Display the first plot
print(plot_group_smrtcell)
# 

# --- Horizontal Line Separator ---
# -------------------------------------------------------


# =======================================================
# --- PLOT 2: Stacked Barplot of 'Gender' within 'SMRTcell_GROUPS' ---
# =======================================================

# 24 bars on the x-axis, 2 colors per bar
plot_gender_smrtcell <- ggplot(metadata, aes(x = SMRTcell_GROUPS, fill = Gender)) +
  geom_bar(position = "stack") +
  
  # Add labels and title
  labs(
    title = "Distribution of 'Gender' within SMRTcell Groups",
    x = "SMRTcell Group",
    y = "Count of Samples",
    fill = "Gender" # Legend title for the fill color
  ) +
  
  # Improve aesthetics: rotate x-axis labels for readability 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  
  # Optional: Use a color-blind friendly palette for two groups 
  scale_fill_manual(values = c("Female" = "darkorchid4", "Male" = "darkorange3"))

# Display the second plot
print(plot_gender_smrtcell)
# 

# --- End of Script ---