
# ============================================================
# TDP-43 STATUS BY NEUROPATHOLOGICAL GROUP
# Final 274-individual PacBio Kinnex cohort
#
# Outputs:
#   - TDP-43 positivity dot plot (PDF and PNG)
#   - TDP-43 status stacked bar plot (PDF and PNG)
#   - Group summary table (CSV)
#   - Sample-level TDP-43 metadata (CSV)
#
# Run directly in Posit Workbench / RStudio
# ============================================================

library(data.table)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# 1. INPUT AND OUTPUT PATHS
# ------------------------------------------------------------

PROJECT <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts"

BASE <- file.path(PROJECT, "isomatch_outputs")

METADATA_FILE <- file.path(
  PROJECT,
  "CWOW_full_long_read_cohort_n276.csv"
)

COUNT_FILE <- file.path(
  BASE,
  "04_confidence_filter_rescuedLR_274",
  "moderate",
  "CWOW_274_rescuedLR.confidence_filtered.counts.tsv.gz"
)

OUT <- file.path(
  BASE,
  "05_DTE_DTU",
  "DTE_DESeq2",
  "moderate",
  "TDP43_plots"
)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. FIGURE SETTINGS
# ------------------------------------------------------------

# Editable figure titles
# No automatic figure numbering

DOT_PLOT_TITLE <- "TDP-43 positivity by neuropathological group"

BAR_PLOT_TITLE <- "TDP-43 status by neuropathological group"

# Small subtitle for the stacked bar plot
BAR_PLOT_SUBTITLE <- paste(
  "Colored = TDP-43 positive; gray = TDP-43 negative.",
  "Missing classifications excluded."
)

# Group colors
group_colors <- c(
  CONTROL = "#377EB8",
  AD_AT   = "#E41A1C",
  LBD_S   = "#6A3D9A",
  LBD_AS  = "#9E9AC8",
  LBD_ATS = "#636363"
)

group_order <- names(group_colors)

# Display labels for manuscript figures
group_labels <- c(
  CONTROL = "CONTROL",
  AD_AT   = "AD(AT)",
  LBD_S   = "LBD(S)",
  LBD_AS  = "LBD(AS)",
  LBD_ATS = "LBD(ATS)"
)

# ------------------------------------------------------------
# 3. READ AND STANDARDIZE METADATA
# ------------------------------------------------------------

meta <- fread(
  METADATA_FILE,
  check.names = FALSE
)

required_columns <- c(
  "Study-Subject-ID",
  "group",
  "TDP.43"
)

missing_columns <- setdiff(
  required_columns,
  names(meta)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Standardize sample identifiers
meta[, sample_id := paste0(
  "CWOW_",
  sub(
    "^CWOW_",
    "",
    trimws(as.character(`Study-Subject-ID`))
  )
)]

# Standardize neuropathological group names
meta[, group := gsub(
  "-",
  "_",
  toupper(trimws(as.character(group)))
)]

# Original metadata uses AD; internal plotting code uses AD_AT
meta[group == "AD", group := "AD_AT"]

cat("\nGroups in original metadata:\n")
print(table(meta$group, useNA = "ifany"))

if (anyDuplicated(meta$sample_id)) {
  stop("Duplicate sample identifiers in metadata.")
}

# ------------------------------------------------------------
# 4. MATCH THE FINAL 274 INDIVIDUALS
# ------------------------------------------------------------

# Read only the count-matrix header
# Avoid loading the large expression matrix

count_header <- names(
  fread(
    COUNT_FILE,
    nrows = 0,
    check.names = FALSE
  )
)

if (
  length(count_header) < 2L ||
  count_header[1] != "transcript_id"
) {
  stop("Expected transcript_id as the first count column.")
}

samples <- count_header[-1L]

if (length(samples) != 274L) {
  stop(
    "Expected 274 individuals; found ",
    length(samples)
  )
}

if (anyDuplicated(samples)) {
  stop("Duplicate sample names in count matrix.")
}

# Match metadata to the exact sample order in the count matrix
meta <- meta[match(samples, sample_id)]

if (anyNA(meta$sample_id)) {
  stop(
    "Missing metadata for samples: ",
    paste(
      samples[is.na(meta$sample_id)],
      collapse = ", "
    )
  )
}

stopifnot(identical(meta$sample_id, samples))

# Validate group labels
unexpected_groups <- setdiff(
  unique(meta$group),
  group_order
)

if (length(unexpected_groups) > 0L) {
  stop(
    "Unexpected group labels: ",
    paste(unexpected_groups, collapse = ", ")
  )
}

meta[, group := factor(
  group,
  levels = group_order
)]

cat("\nFinal cohort group sizes:\n")
print(table(meta$group))

# ------------------------------------------------------------
# 5. STANDARDIZE TDP-43 STATUS
# ------------------------------------------------------------

# Original coding:
#   0 = TDP-43 negative
#   1 = TDP-43 positive
#   NA = missing / unclassified

meta[, TDP43_raw := trimws(
  as.character(`TDP.43`)
)]

meta[
  TDP43_raw %in% c("", "NA", "NaN", "NULL"),
  TDP43_raw := NA_character_
]

unexpected_tdp <- unique(
  meta[
    !is.na(TDP43_raw) &
      !TDP43_raw %in% c("0", "1"),
    TDP43_raw
  ]
)

if (length(unexpected_tdp) > 0L) {
  stop(
    "Unexpected TDP.43 values: ",
    paste(unexpected_tdp, collapse = ", ")
  )
}

meta[, TDP43_status := fifelse(
  TDP43_raw == "1",
  "TDP-43+",
  fifelse(
    TDP43_raw == "0",
    "TDP-43-",
    NA_character_
  )
)]

meta[, TDP43_status := factor(
  TDP43_status,
  levels = c("TDP-43-", "TDP-43+")
)]

cat("\nTDP-43 status by group:\n")
print(table(
  meta$group,
  meta$TDP43_status,
  useNA = "ifany"
))

# ------------------------------------------------------------
# 6. SUMMARIZE TDP-43 DISTRIBUTION
# ------------------------------------------------------------

summary <- meta[, .(
  N_total = .N,
  N_positive = sum(
    TDP43_status == "TDP-43+",
    na.rm = TRUE
  ),
  N_negative = sum(
    TDP43_status == "TDP-43-",
    na.rm = TRUE
  ),
  N_missing = sum(
    is.na(TDP43_status)
  )
), by = group]

summary[, N_assessed := N_positive + N_negative]

summary[, Percent_positive := fifelse(
  N_assessed > 0,
  100 * N_positive / N_assessed,
  NA_real_
)]

summary[, Percent_negative := fifelse(
  N_assessed > 0,
  100 * N_negative / N_assessed,
  NA_real_
)]

summary[, Label := paste0(
  N_positive,
  "/",
  N_assessed
)]

setcolorder(
  summary,
  c(
    "group",
    "N_total",
    "N_assessed",
    "N_positive",
    "N_negative",
    "N_missing",
    "Percent_positive",
    "Percent_negative",
    "Label"
  )
)

setorder(summary, group)

cat("\nTDP-43 summary:\n")
print(summary)

cat("\nTotal individuals:", nrow(meta), "\n")
cat("TDP-43 assessed:", sum(summary$N_assessed), "\n")
cat("TDP-43 missing:", sum(summary$N_missing), "\n")

# Export summary table
fwrite(
  summary,
  file.path(
    OUT,
    "Table_TDP43_by_neuropathology.csv"
  )
)

# Export sample-level classifications
fwrite(
  meta[, .(
    sample_id,
    group,
    TDP43_raw,
    TDP43_status
  )],
  file.path(
    OUT,
    "TDP43_final_cohort_n274.csv"
  )
)

# ------------------------------------------------------------
# 7. CONSISTENT MANUSCRIPT PLOT THEME
# ------------------------------------------------------------

# Figure titles are enabled and editable.
# Axis labels and tick labels use regular-weight font.
# Subtitles are smaller than titles.

theme_manuscript <- function() {
  theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(
        face = "plain",
        color = "black",
        size = 11
      ),

      axis.text.y = element_text(
        face = "plain",
        color = "black",
        size = 11
      ),

      axis.title.x = element_text(
        face = "plain",
        color = "black",
        size = 12
      ),

      axis.title.y = element_text(
        face = "plain",
        color = "black",
        size = 12,
        margin = margin(r = 10)
      ),

      axis.line = element_line(
        color = "black",
        linewidth = 0.5
      ),

      axis.ticks = element_line(
        color = "black",
        linewidth = 0.5
      ),

      plot.title = element_text(
        face = "plain",
        color = "black",
        size = 14,
        hjust = 0,
        margin = margin(b = 5)
      ),

      plot.subtitle = element_text(
        face = "plain",
        color = "#555555",
        size = 9,
        hjust = 0,
        margin = margin(b = 12)
      ),

      legend.text = element_text(
        face = "plain"
      ),

      legend.title = element_text(
        face = "plain"
      ),

      plot.margin = margin(15, 15, 10, 10)
    )
}

# ------------------------------------------------------------
# 8. PLOT 1: TDP-43 POSITIVITY DOT PLOT
# ------------------------------------------------------------

# Color represents neuropathological group.
# Labels represent N positive / N assessed.

# Ensure enough vertical space for count labels
y_max <- max(
  15,
  max(
    summary$Percent_positive,
    na.rm = TRUE
  ) * 1.25
)

p_dot <- ggplot(
  summary,
  aes(
    x = group,
    y = Percent_positive,
    color = group
  )
) +
  geom_segment(
    aes(
      xend = group,
      y = 0,
      yend = Percent_positive
    ),
    linewidth = 0.8,
    alpha = 0.45
  ) +
  geom_point(
    size = 4.5
  ) +
  geom_text(
    aes(label = Label),
    color = "black",
    vjust = -1.2,
    size = 3.8,
    fontface = "plain"
  ) +
  scale_color_manual(
    values = group_colors,
    drop = FALSE
  ) +
  scale_x_discrete(
    labels = group_labels
  ) +
  scale_y_continuous(
    limits = c(0, y_max),
    labels = label_percent(scale = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = DOT_PLOT_TITLE,
    x = NULL,
    y = "TDP-43-positive individuals (%)"
  ) +
  theme_manuscript() +
  theme(
    legend.position = "none"
  )

print(p_dot)

ggsave(
  filename = file.path(
    OUT,
    "TDP43_positive_percentage_by_group.pdf"
  ),
  plot = p_dot,
  width = 7.5,
  height = 5
)

ggsave(
  filename = file.path(
    OUT,
    "TDP43_positive_percentage_by_group.png"
  ),
  plot = p_dot,
  width = 7.5,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 9. PLOT 2: STACKED BAR PLOT OF TDP-43 STATUS
# ------------------------------------------------------------

# TDP-43-positive segments use each group's assigned color.
# TDP-43-negative segments use light gray.
# Missing classifications are excluded from plotted counts.

plot_data <- melt(
  summary,
  id.vars = c(
    "group",
    "N_total",
    "N_assessed",
    "N_missing"
  ),
  measure.vars = c(
    "N_negative",
    "N_positive"
  ),
  variable.name = "Status",
  value.name = "N"
)

plot_data[, Status := factor(
  Status,
  levels = c("N_negative", "N_positive"),
  labels = c("TDP-43-", "TDP-43+")
)]

plot_data[, fill_color := fifelse(
  Status == "TDP-43+",
  unname(group_colors[as.character(group)]),
  "#E5E5E5"
)]

plot_data[, Count_label := fifelse(
  N > 0,
  as.character(N),
  ""
)]

p_bar <- ggplot(
  plot_data,
  aes(
    x = group,
    y = N,
    fill = fill_color,
    group = Status
  )
) +
  geom_col(
    width = 0.7,
    color = "white",
    linewidth = 0.5,
    position = position_stack()
  ) +
  geom_text(
    aes(label = Count_label),
    position = position_stack(vjust = 0.5),
    size = 3.8,
    color = "black",
    fontface = "plain"
  ) +
  scale_fill_identity() +
  scale_x_discrete(
    labels = group_labels
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = BAR_PLOT_TITLE,
    subtitle = BAR_PLOT_SUBTITLE,
    x = NULL,
    y = "Number of individuals"
  ) +
  theme_manuscript()

print(p_bar)

ggsave(
  filename = file.path(
    OUT,
    "TDP43_counts_by_group.pdf"
  ),
  plot = p_bar,
  width = 8,
  height = 5
)

ggsave(
  filename = file.path(
    OUT,
    "TDP43_counts_by_group.png"
  ),
  plot = p_bar,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 10. COMPLETE
# ------------------------------------------------------------

cat("\n====================================\n")
cat("TDP-43 plots completed.\n")
cat("Output directory:\n")
cat(OUT, "\n")
cat("====================================\n")
