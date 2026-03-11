library(RNASeqPower)
library(ggplot2)
library(reshape2)
library(dplyr)
library(PROPER)
library(DESeq2)
library(data.table)

metadata_n580 <- read.delim("/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/metadata/metadata_n580.tsv", header = TRUE, sep = "\t")

rnapower(depth=5, cv=1, effect=c(1.25, 1.5, 1.75, 2),
         alpha= .05, power=c(.8, .9))

#-------------------------------------------------------
# Parameters
sample_sizes <- c(5, 10, 15, 20, 25, 30)  # Sample sizes per group
effect_size <- 2  # Log2 fold change
fdr <- 0.05  # FDR threshold
alpha <- fdr  # Significance level
cv <- 1  # Coefficient of variation (adjust based on your data)

# Calculate power for sequencing depth = 10
power_results <- data.frame(
  Sample_Size = sample_sizes,
  Long_Read = sapply(sample_sizes, function(n) {
    rnapower(n = n, cv = cv, effect = effect_size, alpha = alpha, depth = 5)
  })
)

# Calculate power for sequencing depth = 100
power_results$Short_Read <- sapply(sample_sizes, function(n) {
  rnapower(n = n, cv = cv, effect = effect_size, alpha = alpha, depth = 100)
})

# Reshape data for plotting
power_results_long <- melt(
  power_results, 
  id.vars = "Sample_Size", 
  variable.name = "Coverage_Type", 
  value.name = "Power"
)

# Modify labels for sequencing depth
power_results_long$Coverage_Type <- recode(
  power_results_long$Coverage_Type,
  "Short_Read" = "short-read",
  "Long_Read" = "long-read"
)

# Plot results
ggplot(power_results_long, aes(x = Sample_Size, y = Power, color = Coverage_Type, group = Coverage_Type)) +
  geom_line(aes(linetype = Coverage_Type), size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("short-read" = "blue", "long-read" = "red")) +
  scale_linetype_manual(values = c("short-read" = "solid", "long-read" = "dashed")) +
  labs(
    title = "Power Calculation for RNAseq Isoform Expression",
    x = "Sample Size per Group",
    y = "Power",
    color = "Coverage Type",
    linetype = "Coverage Type"
  ) +
  theme_minimal()
