
# ============================================================
# CWOW LONG-READ DTE: SUMMARY HEATMAPS
#
# Final 274-individual cohort
#
# Generates:
#   1. All 10 pairwise neuropathological comparisons
#   2. Four disease-versus-control comparisons
#
# Counts significantly differentially expressed transcripts:
#   Upregulated:   log2FoldChange > 0
#   Downregulated: log2FoldChange < 0
#   Significance:  FDR-adjusted p < 0.05
#
# Direction is relative to the FIRST group in each comparison.
#
# Uses existing output from:
#   02_DTE_with_TDP43_and_filtering_table.R
#
# Does NOT rerun DESeq2 or modify existing DTE results.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------
# 1. INPUT AND OUTPUT PATHS
# ------------------------------------------------------------

BASE <- paste0(
  "/tgen_labs/jfryer/kolney/LBD_CWOW/",
  "isoform_CWOW/scripts/isomatch_outputs"
)

SCENARIO <- "moderate"

DTE_DIR <- file.path(
  BASE,
  "05_DTE_DTU",
  "DTE_DESeq2",
  SCENARIO
)

OUT <- file.path(
  BASE,
  "05_DTE_DTU",
  "DTE_Figures",
  SCENARIO
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. EDITABLE FIGURE SETTINGS
# ------------------------------------------------------------

# Edit these titles as needed.
# No automatic figure numbering.

ALL_PLOT_TITLE <- paste(
  "Differential transcript expression across",
  "neuropathological comparisons"
)

CONTROL_PLOT_TITLE <- paste(
  "Differential transcript expression in",
  "disease groups versus controls"
)

# Subtitle shared by both figures
PLOT_SUBTITLE <- paste(
  "Significant transcripts (FDR < 0.05).",
  "Direction is relative to the first group",
  "in each comparison."
)

# Significance threshold
FDR <- 0.05

# Optional additional comparison:
# TRUE includes TDP-43+ vs TDP-43- in the full heatmap.
# It is never included in the disease-versus-control heatmap.

INCLUDE_TDP43 <- FALSE

# Heatmap colors
DOWN_COLOR <- "#377EB8"
UP_COLOR   <- "#E41A1C"
ZERO_COLOR <- "#FFFFFF"

# Figure dimensions
ALL_WIDTH <- 9
ALL_HEIGHT <- if (INCLUDE_TDP43) 7 else 6.5

CONTROL_WIDTH <- 9
CONTROL_HEIGHT <- 4.5

# ------------------------------------------------------------
# 3. COMPARISON ORDER AND DISPLAY LABELS
# ------------------------------------------------------------

# These names must match the folders created by 02_DTE.

comparison_labels <- c(
  AD_vs_CONTROL =
    "AD(AT) vs CONTROL",

  LBD_S_vs_CONTROL =
    "LBD(S) vs CONTROL",

  LBD_AS_vs_CONTROL =
    "LBD(AS) vs CONTROL",

  LBD_ATS_vs_CONTROL =
    "LBD(ATS) vs CONTROL",

  LBD_S_vs_AD =
    "LBD(S) vs AD(AT)",

  LBD_AS_vs_AD =
    "LBD(AS) vs AD(AT)",

  LBD_ATS_vs_AD =
    "LBD(ATS) vs AD(AT)",

  LBD_AS_vs_LBD_S =
    "LBD(AS) vs LBD(S)",

  LBD_ATS_vs_LBD_S =
    "LBD(ATS) vs LBD(S)",

  LBD_ATS_vs_LBD_AS =
    "LBD(ATS) vs LBD(AS)"
)

# Four disease-versus-control comparisons
control_comparisons <- c(
  "AD_vs_CONTROL",
  "LBD_S_vs_CONTROL",
  "LBD_AS_vs_CONTROL",
  "LBD_ATS_vs_CONTROL"
)

# Optional TDP-43 comparison
if (INCLUDE_TDP43) {

  comparison_labels <- c(
    comparison_labels,
    TDP43_positive_vs_negative =
      "TDP-43+ vs TDP-43-"
  )

}

# ------------------------------------------------------------
# 4. READ SIGNIFICANT DTE RESULTS
# ------------------------------------------------------------

results_list <- vector(
  "list",
  length(comparison_labels)
)

names(results_list) <- names(comparison_labels)

for (nm in names(comparison_labels)) {

  message("\nReading: ", nm)

  f <- file.path(
    DTE_DIR,
    nm,
    "DTE_FDR05.tsv.gz"
  )

  if (!file.exists(f)) {

    stop(
      "Missing DTE results:\n",
      f,
      "\nCheck whether this comparison completed."
    )

  }

  # Read only the columns needed for summary counts
  d <- fread(
    f,
    select = c(
      "transcript_id",
      "log2FoldChange",
      "padj"
    ),
    showProgress = FALSE
  )

  if (anyDuplicated(d$transcript_id)) {

    stop(
      "Duplicate transcript IDs in:\n",
      f
    )

  }

  # Recheck significance and valid fold changes
  d <- d[
    !is.na(padj) &
      padj < FDR &
      is.finite(log2FoldChange)
  ]

  # Positive log2FC = higher in first group
  n_up <- d[
    log2FoldChange > 0,
    .N
  ]

  # Negative log2FC = lower in first group
  n_down <- d[
    log2FoldChange < 0,
    .N
  ]

  # Significant but exactly zero estimated log2FC
  n_zero <- d[
    log2FoldChange == 0,
    .N
  ]

  results_list[[nm]] <- data.table(
    comparison = nm,
    label = unname(comparison_labels[[nm]]),
    Upregulated = n_up,
    Downregulated = n_down,
    Zero_log2FC = n_zero,
    Total_significant = n_up + n_down + n_zero,
    Total_directional = n_up + n_down
  )

  message(
    "  Upregulated: ", n_up,
    "\n  Downregulated: ", n_down,
    "\n  Zero log2FC: ", n_zero
  )

  rm(d)

}

summary_dt <- rbindlist(
  results_list,
  use.names = TRUE
)

summary_dt[, comparison := factor(
  comparison,
  levels = names(comparison_labels)
)]

setorder(summary_dt, comparison)

cat("\n====================================\n")
cat("DTE SUMMARY\n")
cat("====================================\n")

print(summary_dt)

# ------------------------------------------------------------
# 5. EXPORT SUMMARY TABLES
# ------------------------------------------------------------

# All selected comparisons
fwrite(
  summary_dt,
  file.path(
    OUT,
    "DTE_up_down_summary_all_comparisons.tsv"
  ),
  sep = "\t"
)

# Disease versus control only
control_summary <- copy(
  summary_dt[
    comparison %in% control_comparisons
  ]
)

setorder(control_summary, comparison)

fwrite(
  control_summary,
  file.path(
    OUT,
    "DTE_up_down_summary_disease_vs_control.tsv"
  ),
  sep = "\t"
)

# ------------------------------------------------------------
# 6. PREPARE FULL HEATMAP DATA
# ------------------------------------------------------------

heatmap_dt <- melt(
  summary_dt,
  id.vars = c(
    "comparison",
    "label"
  ),
  measure.vars = c(
    "Upregulated",
    "Downregulated"
  ),
  variable.name = "Direction",
  value.name = "N_DTE"
)

# Column order
heatmap_dt[, Direction := factor(
  Direction,
  levels = c(
    "Upregulated",
    "Downregulated"
  )
)]

# Row order: first comparison appears at top
heatmap_dt[, label := factor(
  label,
  levels = rev(
    unname(comparison_labels)
  )
)]

# Signed values for diverging color scale:
#   Positive = upregulated
#   Negative = downregulated

heatmap_dt[, signed_count := fifelse(
  Direction == "Downregulated",
  -N_DTE,
  N_DTE
)]

# Use the same color scale for both heatmaps
max_count <- max(
  abs(heatmap_dt$signed_count),
  na.rm = TRUE
)

if (max_count == 0) {
  max_count <- 1
}

# Numeric labels
heatmap_dt[, count_label := comma(N_DTE)]

# Text color based on tile intensity
heatmap_dt[, text_color := fifelse(
  abs(signed_count) > 0.55 * max_count,
  "white",
  "black"
)]

# ------------------------------------------------------------
# 7. PREPARE DISEASE-VERSUS-CONTROL DATA
# ------------------------------------------------------------

heatmap_control <- copy(
  heatmap_dt[
    comparison %in% control_comparisons
  ]
)

if (nrow(heatmap_control) != 8L) {

  stop(
    "Expected 8 heatmap cells for four ",
    "disease-versus-control comparisons; found ",
    nrow(heatmap_control)
  )

}

control_labels <- unname(
  comparison_labels[control_comparisons]
)

heatmap_control[, label := factor(
  as.character(label),
  levels = rev(control_labels)
)]

# Both heatmaps use the SAME max_count.

# ------------------------------------------------------------
# 8. CONSISTENT MANUSCRIPT THEME
# ------------------------------------------------------------

theme_manuscript <- function() {

  theme_minimal(base_size = 13) +
    theme(

      # Editable figure titles
      plot.title = element_text(
        face = "plain",
        color = "black",
        size = 14,
        hjust = 0,
        margin = margin(b = 5)
      ),

      # Small subtitle
      plot.subtitle = element_text(
        face = "plain",
        color = "#555555",
        size = 9,
        hjust = 0,
        margin = margin(b = 12)
      ),

      # Regular-weight axis text
      axis.text.x = element_text(
        face = "plain",
        color = "black",
        size = 11
      ),

      axis.text.y = element_text(
        face = "plain",
        color = "black",
        size = 10
      ),

      axis.title.x = element_text(
        face = "plain",
        color = "black",
        size = 12
      ),

      axis.title.y = element_text(
        face = "plain",
        color = "black",
        size = 12
      ),

      axis.ticks = element_blank(),

      panel.grid = element_blank(),

      legend.title = element_text(
        face = "plain",
        color = "black",
        size = 11
      ),

      legend.text = element_text(
        face = "plain",
        color = "black",
        size = 10
      ),

      plot.margin = margin(
        15, 20, 10, 10
      )

    )

}

# ------------------------------------------------------------
# 9. FUNCTION TO GENERATE HEATMAPS
# ------------------------------------------------------------

make_dte_heatmap <- function(
    plot_data,
    plot_title,
    max_count
) {

  ggplot(
    plot_data,
    aes(
      x = Direction,
      y = label,
      fill = signed_count
    )
  ) +

    geom_tile(
      color = "white",
      linewidth = 1,
      width = 0.96,
      height = 0.96
    ) +

    geom_text(
      aes(
        label = count_label,
        color = text_color
      ),
      size = 3.8,
      fontface = "plain"
    ) +

    scale_color_identity() +

    scale_fill_gradient2(
      low = DOWN_COLOR,
      mid = ZERO_COLOR,
      high = UP_COLOR,
      midpoint = 0,
      limits = c(
        -max_count,
        max_count
      ),
      labels = function(x) {
        comma(abs(x))
      },
      name = "Number of DTEs"
    ) +

    scale_x_discrete(
      position = "top",
      drop = FALSE
    ) +

    labs(
      title = plot_title,
      subtitle = PLOT_SUBTITLE,
      x = NULL,
      y = NULL
    ) +

    theme_manuscript() +

    theme(
      legend.position = "right"
    )

}

# ------------------------------------------------------------
# 10. HEATMAP 1: ALL PAIRWISE COMPARISONS
# ------------------------------------------------------------

message("\nGenerating all-comparisons heatmap...")

p_heatmap_all <- make_dte_heatmap(
  plot_data = heatmap_dt,
  plot_title = ALL_PLOT_TITLE,
  max_count = max_count
)

print(p_heatmap_all)

ggsave(
  filename = file.path(
    OUT,
    "DTE_up_down_heatmap_all_comparisons.pdf"
  ),
  plot = p_heatmap_all,
  width = ALL_WIDTH,
  height = ALL_HEIGHT
)

ggsave(
  filename = file.path(
    OUT,
    "DTE_up_down_heatmap_all_comparisons.png"
  ),
  plot = p_heatmap_all,
  width = ALL_WIDTH,
  height = ALL_HEIGHT,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 11. HEATMAP 2: DISEASE VERSUS CONTROL
# ------------------------------------------------------------

message("\nGenerating disease-versus-control heatmap...")

p_heatmap_control <- make_dte_heatmap(
  plot_data = heatmap_control,
  plot_title = CONTROL_PLOT_TITLE,
  max_count = max_count
)

print(p_heatmap_control)

ggsave(
  filename = file.path(
    OUT,
    "DTE_up_down_heatmap_disease_vs_control.pdf"
  ),
  plot = p_heatmap_control,
  width = CONTROL_WIDTH,
  height = CONTROL_HEIGHT
)

ggsave(
  filename = file.path(
    OUT,
    "DTE_up_down_heatmap_disease_vs_control.png"
  ),
  plot = p_heatmap_control,
  width = CONTROL_WIDTH,
  height = CONTROL_HEIGHT,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 12. FINAL VALIDATION
# ------------------------------------------------------------

# Confirm disease-versus-control comparisons
# appear in the expected order.

stopifnot(
  identical(
    as.character(control_summary$comparison),
    control_comparisons
  )
)

# Confirm the two heatmaps use identical color limits.
stopifnot(
  max(abs(heatmap_control$signed_count)) <= max_count
)

# ------------------------------------------------------------
# 13. COMPLETE
# ------------------------------------------------------------

cat("\n====================================\n")
cat("DTE HEATMAPS COMPLETED\n")
cat("====================================\n")

cat("\nGenerated figures:\n")

cat(
  "1. All pairwise comparisons (",
  length(comparison_labels),
  " comparisons)\n",
  sep = ""
)

cat(
  "2. Disease versus CONTROL (4 comparisons)\n"
)

cat("\nOutput directory:\n")
cat(OUT, "\n")

cat("\nSummary tables:\n")
cat("DTE_up_down_summary_all_comparisons.tsv\n")
cat("DTE_up_down_summary_disease_vs_control.tsv\n")

cat("\n====================================\n")
