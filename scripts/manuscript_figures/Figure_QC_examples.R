# ============================================================
# CWOW PACBIO KINNEX
# MANUSCRIPT FIGURES 1 AND 2
#
# 211 CURRENTLY COMPLETED INDIVIDUALS
#
# Disease groups:
#   Control
#   AD
#   LBD-S
#   LBD-AS
#   LBD-ATS
#
# IMPORTANT:
#
# STRUCTURAL-CATEGORY OVERVIEW PANELS use the FULL IsoMatch
# cohort classification so that the following can be shown:
#
#   FSM
#   ISM
#   NIC
#   NNC
#   Fusion
#   Genic
#   Genic intron
#   Antisense
#   Intergenic
#   More junctions
#
# DTU remains based on the separately generated high-confidence
# FSM/NIC/NNC catalog.
#
# FIGURE 1
#
# A. Transcript-model overlap across disease groups (UpSet)
# B. Annotated-gene overlap across disease groups (UpSet)
# C. Detected transcript models per individual
# D. Structural-category composition of cohort transcriptome
# E. Gene/transcript accumulation curves for EACH disease group
#
# FIGURE 2
#
# A. Number of individuals supporting transcripts by category
# B. Number of disease groups supporting transcripts by category
# C. Number of detected transcript models per gene by disease
# D. Transcript length:
#       density by disease
#       stacked histogram by disease/category
#       boxplot by category
# E. Exon number:
#       binned barplot by disease
#       boxplot by disease
#       boxplot by category
# F. Structural-category composition by disease group
#
# BONUS:
# G. Sample transcriptome heatmap/dendrogram
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

library(data.table)
library(Matrix)
library(ggplot2)
library(patchwork)
library(scales)
library(UpSetR)
library(grid)
library(pheatmap)


# ============================================================
# 1. PROJECT PATHS
# ============================================================

outdir <- paste0(
  "/tgen_labs/jfryer/kolney/LBD_CWOW/",
  "isoform_CWOW/isomatch_cohort"
)

figdir <- file.path(
  outdir,
  "manuscript_figures_v2"
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# HC count matrix
# Used for sample-level abundance / dendrogram.
# ------------------------------------------------------------

count_file <- file.path(
  outdir,
  "CWOW_211_high_confidence_ISOMT_sparse_counts.rds"
)


# ------------------------------------------------------------
# Matched 211-sample metadata
# ------------------------------------------------------------

metadata_file <- file.path(
  outdir,
  "CWOW_211_DTU_metadata.tsv"
)


# ------------------------------------------------------------
# FULL IsoMatch classification
#
# 7.77M transcript models
# ------------------------------------------------------------

classification_file <- file.path(
  outdir,
  "CWOW_cohort_gencodev38.classification.txt.gz"
)


# ------------------------------------------------------------
# IsoMatch source PB -> cohort ISOMT mapping
#
# Required to determine which individual expresses which model.
# ------------------------------------------------------------

track_file <- file.path(
  outdir,
  "CWOW_cohort.track.tsv.gz"
)


stopifnot(
  file.exists(count_file),
  file.exists(metadata_file),
  file.exists(classification_file),
  file.exists(track_file)
)


# ============================================================
# 2. GLOBAL PLOT SETTINGS
# ============================================================

# ------------------------------------------------------------
# Disease order
# ------------------------------------------------------------

disease_levels <- c(
  "Control",
  "AD",
  "LBD-S",
  "LBD-AS",
  "LBD-ATS"
)


# ------------------------------------------------------------
# REQUESTED DISEASE COLORS
#
# Control = blue
# AD      = red
# LBD     = purple -> gray
# ------------------------------------------------------------

disease_colors <- c(

  "Control" =
    "#377EB8",

  "AD" =
    "#E41A1C",

  "LBD-S" =
    "#6A3D9A",

  "LBD-AS" =
    "#9E9AC8",

  "LBD-ATS" =
    "#636363"
)


# ------------------------------------------------------------
# Structural-category colors
#
# Chosen to approximately resemble the article.
# ------------------------------------------------------------

category_colors <- c(

  "FSM" =
    "#377EB8",

  "ISM" =
    "#F28E2B",

  "NIC" =
    "#FFC20A",

  "NNC" =
    "#9E9E9E",

  "Fusion" =
    "#6BAED6",

  "Genic" =
    "#70AD47",

  "Genic intron" =
    "#C49A6C",

  "Antisense" =
    "#1F4E79",

  "Intergenic" =
    "#8C2D04",

  "More junctions" =
    "#D45087"
)


category_levels <- names(
  category_colors
)


# ------------------------------------------------------------
# Publication theme
# ------------------------------------------------------------

theme_manuscript <- theme_classic(
  base_size = 12
) +

  theme(

    plot.title =
      element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),

    axis.title =
      element_text(
        face = "bold"
      ),

    axis.text =
      element_text(
        color = "black"
      ),

    legend.title =
      element_text(
        face = "bold"
      ),

    strip.background =
      element_blank(),

    strip.text =
      element_text(
        face = "bold"
      )
  )


# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================


# ============================================================
# 3A. P-VALUE FORMATTER
# ============================================================

format_p <- function(p) {

  if (is.na(p)) {
    return("NA")
  }

  if (p < 2.2e-16) {
    return("P < 2.2 × 10^-16")
  }

  if (p < 0.001) {
    return(
      paste0(
        "P = ",
        format(
          p,
          scientific = TRUE,
          digits = 2
        )
      )
    )
  }

  paste0(
    "P = ",
    signif(
      p,
      3
    )
  )
}


# ============================================================
# 3B. ALL PAIRWISE WILCOXON TESTS
# ============================================================

pairwise_wilcox_table <- function(
  data,
  value_col,
  group_col
) {

  groups <- disease_levels

  pairs <- combn(
    groups,
    2,
    simplify = FALSE
  )

  results <- lapply(
    pairs,
    function(z) {

      x <- data[
        get(group_col) == z[1],
        get(value_col)
      ]

      y <- data[
        get(group_col) == z[2],
        get(value_col)
      ]

      test <- tryCatch(
        wilcox.test(
          x,
          y,
          exact = FALSE
        ),
        error = function(e) NULL
      )

      data.table(
        group1 = z[1],
        group2 = z[2],
        p_value =
          ifelse(
            is.null(test),
            NA_real_,
            test$p.value
          )
      )
    }
  )

  out <- rbindlist(
    results
  )

  out[
    ,
    FDR :=
      p.adjust(
        p_value,
        method = "BH"
      )
  ]

  out[
    ,
    label :=
      vapply(
        FDR,
        format_p,
        character(1)
      )
  ]

  out
}


# ============================================================
# 3C. ADD CONTROL COMPARISON BRACKETS
#
# To avoid unreadable figures, all 10 pairwise comparisons
# are saved to TSV but only disease-vs-Control comparisons
# are drawn on the plot.
# ============================================================

add_control_brackets <- function(
  p,
  stats,
  ymax,
  increment
) {

  control_stats <- stats[
    group1 == "Control" |
      group2 == "Control"
  ]

  control_stats <- control_stats[
    group1 != group2
  ]

  other_group <- ifelse(
    control_stats$group1 == "Control",
    control_stats$group2,
    control_stats$group1
  )

  control_stats[
    ,
    other :=
      other_group
  ]

  control_stats[
    ,
    x1 :=
      match(
        "Control",
        disease_levels
      )
  ]

  control_stats[
    ,
    x2 :=
      match(
        other,
        disease_levels
      )
  ]

  control_stats[
    ,
    ypos :=
      ymax +
      seq_len(.N) *
      increment
  ]


  for (i in seq_len(nrow(control_stats))) {

    row <- control_stats[i]

    p <- p +

      annotate(
        "segment",
        x = row$x1,
        xend = row$x2,
        y = row$ypos,
        yend = row$ypos,
        linewidth = 0.4
      ) +

      annotate(
        "segment",
        x = row$x1,
        xend = row$x1,
        y = row$ypos,
        yend = row$ypos - increment * 0.15,
        linewidth = 0.4
      ) +

      annotate(
        "segment",
        x = row$x2,
        xend = row$x2,
        y = row$ypos,
        yend = row$ypos - increment * 0.15,
        linewidth = 0.4
      ) +

      annotate(
        "text",
        x = mean(
          c(
            row$x1,
            row$x2
          )
        ),
        y = row$ypos +
          increment * 0.07,
        label = row$label,
        size = 3
      )
  }

  p
}


# ============================================================
# 4. LOAD HC COUNT MATRIX
# ============================================================

cat(
  "\nLoading HC count matrix...\n"
)

counts <- readRDS(
  count_file
)

cat(
  "HC matrix: ",
  format(
    nrow(counts),
    big.mark = ","
  ),
  " transcripts x ",
  ncol(counts),
  " individuals\n",
  sep = ""
)


# ============================================================
# 5. LOAD METADATA
# ============================================================

cat(
  "\nLoading metadata...\n"
)

meta <- fread(
  metadata_file
)


# ------------------------------------------------------------
# Make robust individual ID
# ------------------------------------------------------------

if (
  !"individual_id" %in%
    names(meta)
) {

  meta[
    ,
    individual_id :=
      sub(
        "-.*$",
        "",
        sub(
          "^CWOW_",
          "",
          sample
        )
      )
  ]
}


meta[
  ,
  disease_plot :=
    fcase(

      as.character(disease) ==
        "CONTROL",
      "Control",

      as.character(disease) ==
        "AD",
      "AD",

      as.character(disease) ==
        "LBD_S",
      "LBD-S",

      as.character(disease) ==
        "LBD_AS",
      "LBD-AS",

      as.character(disease) ==
        "LBD_ATS",
      "LBD-ATS",

      default =
        NA_character_
    )
]


meta[
  ,
  disease_plot :=
    factor(
      disease_plot,
      levels =
        disease_levels
    )
]


cat(
  "\nDisease sample sizes:\n"
)

print(
  meta[
    ,
    .N,
    by = disease_plot
  ]
)


# ============================================================
# 6. LOAD FULL ISOMATCH CLASSIFICATION
# ============================================================

cat(
  "\nReading full IsoMatch classification...\n"
)


classification <- fread(
  cmd = paste(
    "zcat",
    shQuote(
      classification_file
    )
  ),
  select = c(
    "isoform_id",
    "query_length",
    "query_exon_n",
    "structural_category",
    "subcategory",
    "ref_gene_id",
    "ref_gene_name"
  )
)


setnames(
  classification,
  "isoform_id",
  "transcript_id"
)


cat(
  "Full cohort transcript models: ",
  format(
    nrow(classification),
    big.mark = ","
  ),
  "\n",
  sep = ""
)


# ============================================================
# 7. STRUCTURAL CATEGORY LABELS
# ============================================================

classification[
  ,
  category :=
    fcase(

      structural_category ==
        "full-splice_match",
      "FSM",

      structural_category ==
        "incomplete-splice_match",
      "ISM",

      structural_category ==
        "novel_in_catalog",
      "NIC",

      structural_category ==
        "novel_not_in_catalog",
      "NNC",

      structural_category ==
        "fusion",
      "Fusion",

      structural_category ==
        "genic",
      "Genic",

      structural_category ==
        "genic_intron",
      "Genic intron",

      structural_category ==
        "antisense",
      "Antisense",

      structural_category ==
        "intergenic",
      "Intergenic",

      structural_category ==
        "moreJunctions",
      "More junctions",

      default =
        structural_category
    )
]


classification[
  ,
  category :=
    factor(
      category,
      levels =
        category_levels
    )
]


cat(
  "\nFull structural-category summary:\n"
)

print(
  classification[
    ,
    .N,
    by = category
  ][
    order(
      category
    )
  ]
)


# ============================================================
# 8. LOAD ISOMATCH TRACK FILE
#
# This tells us which individual contributed each cohort model.
# ============================================================

cat(
  "\nReading IsoMatch track file...\n"
)


track <- fread(
  cmd = paste(
    "zcat",
    shQuote(
      track_file
    )
  ),
  select = c(
    "merged_tx_id",
    "src_file_name"
  )
)


setnames(
  track,
  "merged_tx_id",
  "transcript_id"
)


# ------------------------------------------------------------
# Stable individual ID
# ------------------------------------------------------------

track[
  ,
  individual_id :=
    sub(
      "-.*$",
      "",
      sub(
        "^CWOW_",
        "",
        sub(
          "\\.filtered\\.gtf$",
          "",
          src_file_name
        )
      )
    )
]


# ------------------------------------------------------------
# One transcript can have multiple source PB models within
# the same individual. Collapse to transcript x individual.
# ------------------------------------------------------------

presence <- unique(
  track[
    ,
    .(
      transcript_id,
      individual_id
    )
  ]
)


rm(track)

gc()


# ============================================================
# 9. ADD DISEASE TO PRESENCE TABLE
# ============================================================

meta_map <- unique(
  meta[
    ,
    .(
      individual_id,
      disease_plot
    )
  ]
)


presence <- merge(
  presence,
  meta_map,
  by = "individual_id",
  all.x = TRUE
)


presence <- presence[
  !is.na(
    disease_plot
  )
]


cat(
  "\nTranscript-individual observations:",
  format(
    nrow(presence),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Individuals represented:",
  uniqueN(
    presence$individual_id
  ),
  "\n"
)


# ============================================================
# 10. TRANSCRIPT x DISEASE TABLE
# ============================================================

tx_group <- unique(
  presence[
    ,
    .(
      transcript_id,
      disease_plot
    )
  ]
)


# ============================================================
# 11. GENE ANNOTATION LOOKUP
# ============================================================

gene_lookup <- classification[
  !is.na(ref_gene_id) &
    ref_gene_id != "",
  .(
    transcript_id,
    gene_id =
      ref_gene_id,
    gene_name =
      ref_gene_name
  )
]


gene_lookup <- unique(
  gene_lookup,
  by = "transcript_id"
)


# ============================================================
# 12. GENE x DISEASE TABLE
# ============================================================

gene_group <- merge(
  tx_group,
  gene_lookup,
  by = "transcript_id",
  all = FALSE
)


gene_group <- unique(
  gene_group[
    ,
    .(
      gene_id,
      disease_plot
    )
  ]
)


# ============================================================
# FIGURE 1A
# TRANSCRIPT UPSET
# ============================================================

cat(
  "\nCreating transcript intersection combinations...\n"
)


tx_combo <- tx_group[
  ,
  .(
    combination =
      paste(
        sort(
          as.character(
            disease_plot
          )
        ),
        collapse = "&"
      )
  ),
  by = transcript_id
][
  ,
  .N,
  by = combination
]


# ------------------------------------------------------------
# UpSetR fromExpression works best with syntactically simple
# names. Temporarily replace hyphens with underscores.
# ------------------------------------------------------------

tx_combo[
  ,
  combo_safe :=
    gsub(
      "LBD-S",
      "LBD_S",
      combination,
      fixed = TRUE
    )
]

tx_combo[
  ,
  combo_safe :=
    gsub(
      "LBD-AS",
      "LBD_AS",
      combo_safe,
      fixed = TRUE
    )
]

tx_combo[
  ,
  combo_safe :=
    gsub(
      "LBD-ATS",
      "LBD_ATS",
      combo_safe,
      fixed = TRUE
    )
]


tx_upset_data <- UpSetR::fromExpression(
  setNames(
    tx_combo$N,
    tx_combo$combo_safe
  )
)


names(tx_upset_data)[
  names(tx_upset_data) ==
    "LBD_S"
] <- "LBD-S"

names(tx_upset_data)[
  names(tx_upset_data) ==
    "LBD_AS"
] <- "LBD-AS"

names(tx_upset_data)[
  names(tx_upset_data) ==
    "LBD_ATS"
] <- "LBD-ATS"


# ------------------------------------------------------------
# Save PNG
# ------------------------------------------------------------

png(
  filename =
    file.path(
      figdir,
      "Figure1A_transcript_UpSet.png"
    ),
  width = 10,
  height = 6,
  units = "in",
  res = 600
)


UpSetR::upset(

  tx_upset_data,

  sets =
    disease_levels,

  keep.order =
    TRUE,

  order.by =
    "freq",

  nintersects =
    NA,

  mb.ratio =
    c(
      0.65,
      0.35
    ),

  mainbar.y.label =
    "Transcript models",

  sets.x.label =
    "Transcript models per group",

  text.scale =
    c(
      1.3,
      1.2,
      1.0,
      1.0,
      1.1,
      1.0
    )
)


dev.off()


# ------------------------------------------------------------
# PDF
# ------------------------------------------------------------

pdf(
  file =
    file.path(
      figdir,
      "Figure1A_transcript_UpSet.pdf"
    ),
  width = 10,
  height = 6
)


UpSetR::upset(

  tx_upset_data,

  sets =
    disease_levels,

  keep.order =
    TRUE,

  order.by =
    "freq",

  nintersects =
    NA,

  mainbar.y.label =
    "Transcript models",

  sets.x.label =
    "Transcript models per group"
)


dev.off()


# ============================================================
# FIGURE 1B
# GENE UPSET
# ============================================================

gene_combo <- gene_group[
  ,
  .(
    combination =
      paste(
        sort(
          as.character(
            disease_plot
          )
        ),
        collapse = "&"
      )
  ),
  by = gene_id
][
  ,
  .N,
  by = combination
]


gene_combo[
  ,
  combo_safe :=
    gsub(
      "LBD-S",
      "LBD_S",
      combination,
      fixed = TRUE
    )
]

gene_combo[
  ,
  combo_safe :=
    gsub(
      "LBD-AS",
      "LBD_AS",
      combo_safe,
      fixed = TRUE
    )
]

gene_combo[
  ,
  combo_safe :=
    gsub(
      "LBD-ATS",
      "LBD_ATS",
      combo_safe,
      fixed = TRUE
    )
]


gene_upset_data <- UpSetR::fromExpression(
  setNames(
    gene_combo$N,
    gene_combo$combo_safe
  )
)


names(gene_upset_data)[
  names(gene_upset_data) ==
    "LBD_S"
] <- "LBD-S"

names(gene_upset_data)[
  names(gene_upset_data) ==
    "LBD_AS"
] <- "LBD-AS"

names(gene_upset_data)[
  names(gene_upset_data) ==
    "LBD_ATS"
] <- "LBD-ATS"


png(
  filename =
    file.path(
      figdir,
      "Figure1B_gene_UpSet.png"
    ),
  width = 10,
  height = 6,
  units = "in",
  res = 600
)


UpSetR::upset(

  gene_upset_data,

  sets =
    disease_levels,

  keep.order =
    TRUE,

  order.by =
    "freq",

  nintersects =
    NA,

  mainbar.y.label =
    "Annotated genes",

  sets.x.label =
    "Genes per group"
)


dev.off()


pdf(
  file =
    file.path(
      figdir,
      "Figure1B_gene_UpSet.pdf"
    ),
  width = 10,
  height = 6
)


UpSetR::upset(

  gene_upset_data,

  sets =
    disease_levels,

  keep.order =
    TRUE,

  order.by =
    "freq",

  nintersects =
    NA,

  mainbar.y.label =
    "Annotated genes",

  sets.x.label =
    "Genes per group"
)


dev.off()


# ============================================================
# FIGURE 1C
# NUMBER OF TRANSCRIPT MODELS DETECTED PER INDIVIDUAL
# ============================================================

sample_tx_counts <- presence[
  ,
  .(
    N_transcripts =
      uniqueN(
        transcript_id
      )
  ),
  by = .(
    individual_id,
    disease_plot
  )
]


sample_tx_counts[
  ,
  disease_plot :=
    factor(
      disease_plot,
      levels =
        disease_levels
    )
]


p1C_stats <- pairwise_wilcox_table(
  sample_tx_counts,
  "N_transcripts",
  "disease_plot"
)


fwrite(
  p1C_stats,
  file.path(
    figdir,
    "Figure1C_pairwise_Wilcoxon.tsv"
  ),
  sep = "\t"
)


p1C <- ggplot(
  sample_tx_counts,
  aes(
    x = disease_plot,
    y = N_transcripts,
    fill = disease_plot
  )
) +

  geom_boxplot(
    width = 0.65,
    outlier.shape = NA
  ) +

  geom_jitter(
    aes(
      color = disease_plot
    ),
    width = 0.15,
    size = 1.5,
    alpha = 0.55
  ) +

  scale_fill_manual(
    values =
      disease_colors
  ) +

  scale_color_manual(
    values =
      disease_colors
  ) +

  scale_y_continuous(
    labels =
      comma,
    expand =
      expansion(
        mult =
          c(
            0.03,
            0.25
          )
      )
  ) +

  labs(
    title =
      "C",
    x =
      NULL,
    y =
      "Number of detected transcript models"
  ) +

  guides(
    fill = "none",
    color = "none"
  ) +

  theme_manuscript


p1C <- add_control_brackets(

  p =
    p1C,

  stats =
    p1C_stats,

  ymax =
    max(
      sample_tx_counts$N_transcripts
    ),

  increment =
    max(
      sample_tx_counts$N_transcripts
    ) * 0.045
)


# ============================================================
# FIGURE 1D
# FULL STRUCTURAL CATEGORY COMPOSITION
# ============================================================

category_summary <- classification[
  ,
  .(
    N =
      .N
  ),
  by = category
]


category_summary[
  ,
  percent :=
    100 *
    N /
    sum(N)
]


category_summary[
  ,
  label :=
    paste0(
      category,
      "\n",
      round(
        percent,
        1
      ),
      "%"
    )
]


p1D <- ggplot(
  category_summary,
  aes(
    x = "",
    y = N,
    fill = category
  )
) +

  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.2
  ) +

  coord_polar(
    theta = "y"
  ) +

  scale_fill_manual(
    values =
      category_colors,
    drop = FALSE
  ) +

  labs(
    title =
      "D",
    fill =
      NULL
  ) +

  theme_void(
    base_size = 11
  ) +

  theme(
    plot.title =
      element_text(
        face = "bold",
        hjust = 0
      )
  )


# ============================================================
# FIGURE 1E
# INDIVIDUAL-ACCUMULATION / RAREFACTION CURVES
#
# NOTE:
# This is sample accumulation, not read-level subsampling.
#
# For final manuscript figures you may increase N_PERM from
# 5 to 20.
# ============================================================

cat(
  "\nCalculating disease-specific accumulation curves...\n"
)


set.seed(
  12345
)


N_PERM <- 5


classification_gene <- classification[
  ,
  .(
    transcript_id,
    gene_id =
      ref_gene_id
  )
]


setkey(
  classification_gene,
  transcript_id
)


rarefaction_one_group <- function(
  group_value,
  n_perm = 5
) {

  group_samples <- meta[
    disease_plot == group_value,
    individual_id
  ]


  gp <- presence[
    disease_plot == group_value,
    .(
      transcript_id,
      individual_id
    )
  ]


  gp <- classification_gene[
    gp,
    on = "transcript_id"
  ]


  output <- vector(
    "list",
    n_perm
  )


  for (p in seq_len(n_perm)) {

    permuted <- sample(
      group_samples
    )


    rank_table <- data.table(

      individual_id =
        permuted,

      rank =
        seq_along(
          permuted
        )
    )


    x <- rank_table[
      gp,
      on = "individual_id",
      nomatch = 0
    ]


    first_tx <- x[
      ,
      .(
        first_seen =
          min(rank)
      ),
      by = transcript_id
    ]


    tx_new <- tabulate(
      first_tx$first_seen,
      nbins =
        length(
          group_samples
        )
    )


    first_gene <- x[
      !is.na(gene_id) &
        gene_id != "",
      .(
        first_seen =
          min(rank)
      ),
      by = gene_id
    ]


    gene_new <- tabulate(
      first_gene$first_seen,
      nbins =
        length(
          group_samples
        )
    )


    output[[p]] <- data.table(

      disease =
        group_value,

      permutation =
        p,

      N_samples =
        seq_along(
          group_samples
        ),

      N_transcripts =
        cumsum(
          tx_new
        ),

      N_genes =
        cumsum(
          gene_new
        )
    )
  }


  rbindlist(
    output
  )
}


rarefaction_raw <- rbindlist(

  lapply(
    disease_levels,
    rarefaction_one_group,
    n_perm =
      N_PERM
  )
)


rarefaction_summary <- rarefaction_raw[
  ,
  .(
    mean_transcripts =
      mean(
        N_transcripts
      ),

    mean_genes =
      mean(
        N_genes
      )
  ),
  by = .(
    disease,
    N_samples
  )
]


rarefaction_long <- rbindlist(
  list(

    rarefaction_summary[
      ,
      .(
        disease,
        N_samples,
        Type =
          "Transcript",
        N =
          mean_transcripts
      )
    ],

    rarefaction_summary[
      ,
      .(
        disease,
        N_samples,
        Type =
          "Gene",
        N =
          mean_genes
      )
    ]
  )
)


rarefaction_long[
  ,
  disease :=
    factor(
      disease,
      levels =
        disease_levels
    )
]


p1E <- ggplot(
  rarefaction_long,
  aes(
    x = N_samples,
    y = N,
    color = disease,
    linetype = Type
  )
) +

  geom_line(
    linewidth = 0.9
  ) +

  facet_wrap(
    ~ disease,
    scales = "free_x",
    nrow = 2
  ) +

  scale_color_manual(
    values =
      disease_colors
  ) +

  scale_linetype_manual(
    values = c(
      "Gene" = "solid",
      "Transcript" = "dashed"
    )
  ) +

  scale_y_continuous(
    labels =
      comma
  ) +

  labs(
    title =
      "E",
    x =
      "Number of individuals",
    y =
      "Cumulative detected genes/transcripts",
    color =
      NULL,
    linetype =
      NULL
  ) +

  theme_manuscript +

  theme(
    legend.position =
      "bottom"
  )


# ============================================================
# SAVE FIGURE 1 C-E
# ============================================================

ggsave(
  file.path(
    figdir,
    "Figure1C_transcripts_per_individual.png"
  ),
  p1C,
  width = 6.5,
  height = 5,
  dpi = 600
)


ggsave(
  file.path(
    figdir,
    "Figure1D_structural_categories.png"
  ),
  p1D,
  width = 6.5,
  height = 5,
  dpi = 600
)


ggsave(
  file.path(
    figdir,
    "Figure1E_disease_accumulation.png"
  ),
  p1E,
  width = 13,
  height = 8,
  dpi = 600
)


# ============================================================
# FIGURE 2A
# NUMBER OF INDIVIDUALS SUPPORTING EACH TRANSCRIPT
# ============================================================

cat(
  "\nCalculating transcript recurrence...\n"
)


tx_support <- presence[
  ,
  .(
    N_samples =
      uniqueN(
        individual_id
      )
  ),
  by = transcript_id
]


tx_support <- classification[
  tx_support,
  on = "transcript_id",
  nomatch = 0
]


tx_support[
  ,
  support_bin :=
    cut(

      N_samples,

      breaks = c(
        0,
        1,
        4,
        9,
        19,
        49,
        99,
        Inf
      ),

      labels = c(
        "1",
        "2–4",
        "5–9",
        "10–19",
        "20–49",
        "50–99",
        "≥100"
      )
    )
]


support_summary <- tx_support[
  ,
  .N,
  by = .(
    category,
    support_bin
  )
]


support_summary[
  ,
  proportion :=
    N /
    sum(N),
  by = category
]


p2A <- ggplot(
  support_summary,
  aes(
    x = category,
    y = proportion,
    fill = support_bin
  )
) +

  geom_col(
    width = 0.82
  ) +

  scale_y_continuous(
    labels =
      percent_format(
        accuracy = 1
      )
  ) +

  labs(
    title =
      "A",
    x =
      NULL,
    y =
      "Count (%)",
    fill =
      "Number of individuals\nexpressing transcript"
  ) +

  theme_manuscript +

  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


# ============================================================
# FIGURE 2B
# NUMBER OF DISEASE GROUPS IN WHICH EACH TRANSCRIPT OCCURS
# ============================================================

tx_disease_support <- tx_group[
  ,
  .(
    N_disease_groups =
      uniqueN(
        disease_plot
      )
  ),
  by = transcript_id
]


tx_disease_support <- classification[
  tx_disease_support,
  on = "transcript_id",
  nomatch = 0
]


group_support_summary <- tx_disease_support[
  ,
  .N,
  by = .(
    category,
    N_disease_groups
  )
]


group_support_summary[
  ,
  proportion :=
    N /
    sum(N),
  by = category
]


p2B <- ggplot(
  group_support_summary,
  aes(
    x = category,
    y = proportion,
    fill = factor(
      N_disease_groups
    )
  )
) +

  geom_col(
    width = 0.82
  ) +

  scale_y_continuous(
    labels =
      percent_format(
        accuracy = 1
      )
  ) +

  labs(
    title =
      "B",
    x =
      NULL,
    y =
      "Count (%)",
    fill =
      "Number of disease groups\nexpressing transcript"
  ) +

  theme_manuscript +

  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


# ============================================================
# PREPARE TRANSCRIPT x GROUP ANNOTATION FOR C-F
# ============================================================

cat(
  "\nJoining disease presence to transcript annotation...\n"
)


tx_group_annotation <- classification[
  tx_group,
  on = "transcript_id",
  nomatch = 0
]


tx_group_annotation[
  ,
  disease_plot :=
    factor(
      disease_plot,
      levels =
        disease_levels
    )
]


# ============================================================
# FIGURE 2C
# NUMBER OF ISOFORMS PER GENE BY DISEASE
# ============================================================

gene_tx_by_group <- tx_group_annotation[
  !is.na(ref_gene_id) &
    ref_gene_id != "",
  .(
    N_isoforms =
      uniqueN(
        transcript_id
      )
  ),
  by = .(
    disease_plot,
    ref_gene_id
  )
]


p2C_stats <- pairwise_wilcox_table(
  gene_tx_by_group,
  "N_isoforms",
  "disease_plot"
)


fwrite(
  p2C_stats,
  file.path(
    figdir,
    "Figure2C_pairwise_Wilcoxon.tsv"
  ),
  sep = "\t"
)


p2C <- ggplot(
  gene_tx_by_group,
  aes(
    x = disease_plot,
    y = N_isoforms,
    fill = disease_plot
  )
) +

  geom_boxplot(
    width = 0.65,
    outlier.alpha = 0.20,
    outlier.size = 0.5
  ) +

  scale_fill_manual(
    values =
      disease_colors
  ) +

  scale_y_continuous(
    trans =
      "log10",
    expand =
      expansion(
        mult =
          c(
            0.02,
            0.25
          )
      )
  ) +

  labs(
    title =
      "C",
    x =
      NULL,
    y =
      "Number of transcript models per gene"
  ) +

  guides(
    fill =
      "none"
  ) +

  theme_manuscript


p2C <- add_control_brackets(

  p =
    p2C,

  stats =
    p2C_stats,

  ymax =
    max(
      gene_tx_by_group$N_isoforms
    ),

  increment =
    max(
      gene_tx_by_group$N_isoforms
    ) * 0.12
)


# ============================================================
# FIGURE 2D
# TRANSCRIPT LENGTH
# ============================================================


# ------------------------------------------------------------
# D1. Density by disease
# ------------------------------------------------------------

length_99 <- quantile(
  tx_group_annotation$query_length,
  0.99,
  na.rm = TRUE
)


p2D_density <- ggplot(
  tx_group_annotation[
    query_length <=
      length_99
  ],
  aes(
    x = query_length,
    color = disease_plot
  )
) +

  geom_density(
    linewidth = 0.9
  ) +

  scale_color_manual(
    values =
      disease_colors
  ) +

  labs(
    title =
      "D",
    x =
      "Transcript length (bp)",
    y =
      "Density",
    color =
      NULL
  ) +

  theme_manuscript


# ------------------------------------------------------------
# D2. Stacked histograms by disease/category
# ------------------------------------------------------------

p2D_hist <- ggplot(
  tx_group_annotation[
    query_length <=
      length_99
  ],
  aes(
    x = query_length,
    fill = category
  )
) +

  geom_histogram(
    bins = 40
  ) +

  facet_wrap(
    ~ disease_plot,
    nrow = 1,
    scales = "free_y"
  ) +

  scale_fill_manual(
    values =
      category_colors,
    drop = FALSE
  ) +

  labs(
    x =
      "Transcript length (bp)",
    y =
      "Count",
    fill =
      "Category"
  ) +

  theme_manuscript +

  theme(
    legend.position =
      "bottom"
  )


# ------------------------------------------------------------
# D3. Length boxplot by category
# ------------------------------------------------------------

p2D_category <- ggplot(
  classification[
    query_length <=
      quantile(
        query_length,
        0.995,
        na.rm = TRUE
      )
  ],
  aes(
    x = category,
    y = query_length,
    fill = category
  )
) +

  geom_boxplot(
    outlier.alpha = 0.15,
    outlier.size = 0.4
  ) +

  scale_fill_manual(
    values =
      category_colors,
    drop = FALSE
  ) +

  labs(
    x =
      NULL,
    y =
      "Transcript length (bp)"
  ) +

  guides(
    fill =
      "none"
  ) +

  theme_manuscript +

  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


# ------------------------------------------------------------
# Assemble panel D
# ------------------------------------------------------------

p2D <- (
  p2D_density |
    p2D_category
) /
  p2D_hist +
  plot_layout(
    heights =
      c(
        1,
        1.15
      )
  )


# ============================================================
# FIGURE 2E
# EXON NUMBER
# ============================================================


# ------------------------------------------------------------
# Exon-number bins
# ------------------------------------------------------------

tx_group_annotation[
  ,
  exon_bin :=
    cut(
      query_exon_n,
      breaks = c(
        0,
        1,
        4,
        9,
        19,
        Inf
      ),
      labels = c(
        "1",
        "2–4",
        "5–9",
        "10–19",
        "≥20"
      )
    )
]


exon_bar <- tx_group_annotation[
  ,
  .N,
  by = .(
    disease_plot,
    exon_bin
  )
]


exon_bar[
  ,
  proportion :=
    N /
    sum(N),
  by = disease_plot
]


# ------------------------------------------------------------
# E1. Grouped percentage barplot
# ------------------------------------------------------------

p2E_bar <- ggplot(
  exon_bar,
  aes(
    x = exon_bin,
    y = proportion,
    fill = disease_plot
  )
) +

  geom_col(
    position =
      position_dodge(
        width = 0.85
      ),
    width = 0.8
  ) +

  scale_fill_manual(
    values =
      disease_colors
  ) +

  scale_y_continuous(
    labels =
      percent_format(
        accuracy = 1
      )
  ) +

  labs(
    title =
      "E",
    x =
      "Number of exons per transcript",
    y =
      "Count (%)",
    fill =
      NULL
  ) +

  theme_manuscript


# ------------------------------------------------------------
# E2. Boxplot by disease
# ------------------------------------------------------------

p2E_disease <- ggplot(
  tx_group_annotation[
    query_exon_n <=
      quantile(
        query_exon_n,
        0.995,
        na.rm = TRUE
      )
  ],
  aes(
    x = disease_plot,
    y = query_exon_n,
    fill = disease_plot
  )
) +

  geom_boxplot(
    outlier.alpha = 0.15,
    outlier.size = 0.4
  ) +

  scale_fill_manual(
    values =
      disease_colors
  ) +

  labs(
    x =
      NULL,
    y =
      "Number of exons per transcript"
  ) +

  guides(
    fill =
      "none"
  ) +

  theme_manuscript


# ------------------------------------------------------------
# E3. Boxplot by structural category
# ------------------------------------------------------------

p2E_category <- ggplot(
  classification[
    query_exon_n <=
      quantile(
        query_exon_n,
        0.995,
        na.rm = TRUE
      )
  ],
  aes(
    x = category,
    y = query_exon_n,
    fill = category
  )
) +

  geom_boxplot(
    outlier.alpha = 0.15,
    outlier.size = 0.4
  ) +

  scale_fill_manual(
    values =
      category_colors,
    drop = FALSE
  ) +

  labs(
    x =
      NULL,
    y =
      "Number of exons per transcript"
  ) +

  guides(
    fill =
      "none"
  ) +

  theme_manuscript +

  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


p2E <- (
  p2E_bar |
    p2E_disease |
    p2E_category
)


# ============================================================
# FIGURE 2F
# CATEGORY COMPOSITION BY DISEASE
# ============================================================

category_by_group <- tx_group_annotation[
  ,
  .(
    N =
      uniqueN(
        transcript_id
      )
  ),
  by = .(
    disease_plot,
    category
  )
]


category_by_group[
  ,
  proportion :=
    N /
    sum(N),
  by = disease_plot
]


p2F <- ggplot(
  category_by_group,
  aes(
    x = disease_plot,
    y = proportion,
    fill = category
  )
) +

  geom_col(
    width = 0.75
  ) +

  scale_fill_manual(
    values =
      category_colors,
    drop = FALSE
  ) +

  scale_y_continuous(
    labels =
      percent_format(
        accuracy = 1
      )
  ) +

  labs(
    title =
      "F",
    x =
      NULL,
    y =
      "Count (%)",
    fill =
      "Structural category"
  ) +

  theme_manuscript +

  theme(
    legend.position =
      "right"
  )


# ============================================================
# FIGURE 2G
# SAMPLE TRANSCRIPTOME HEATMAP + DENDROGRAM
#
# Uses HC count matrix.
#
# Select 2,000 abundant transcripts, then 500 most variable.
# ============================================================

cat(
  "\nCreating sample clustering heatmap...\n"
)


library_sizes <- Matrix::colSums(
  counts
)


# ------------------------------------------------------------
# Total transcript abundance
# ------------------------------------------------------------

tx_total <- Matrix::rowSums(
  counts
)


top_abundant <- order(
  tx_total,
  decreasing = TRUE
)[
  seq_len(
    min(
      2000,
      length(
        tx_total
      )
    )
  )
]


counts_top <- counts[
  top_abundant,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# CPM for only 2,000 transcripts
# ------------------------------------------------------------

cpm_top <- sweep(
  as.matrix(
    counts_top
  ),
  2,
  library_sizes,
  "/"
) *
  1e6


log_cpm_top <- log2(
  cpm_top +
    1
)


# ------------------------------------------------------------
# Select most variable 500
# ------------------------------------------------------------

tx_var <- apply(
  log_cpm_top,
  1,
  var
)


top_variable <- order(
  tx_var,
  decreasing = TRUE
)[
  seq_len(
    min(
      500,
      length(
        tx_var
      )
    )
  )
]


heatmap_matrix <- log_cpm_top[
  top_variable,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Row z-score
# ------------------------------------------------------------

heatmap_matrix <- t(
  scale(
    t(
      heatmap_matrix
    )
  )
)


heatmap_matrix[
  is.na(
    heatmap_matrix
  )
] <- 0


# ------------------------------------------------------------
# Heatmap metadata
# ------------------------------------------------------------

heatmap_meta <- meta[
  match(
    colnames(
      heatmap_matrix
    ),
    sample
  ),
  .(
    Disease =
      disease_plot
  )
]


heatmap_meta <- as.data.frame(
  heatmap_meta
)


rownames(
  heatmap_meta
) <- colnames(
  heatmap_matrix
)


annotation_colors <- list(
  Disease =
    disease_colors
)


pheatmap(

  heatmap_matrix,

  cluster_rows =
    TRUE,

  cluster_cols =
    TRUE,

  show_rownames =
    FALSE,

  show_colnames =
    FALSE,

  annotation_col =
    heatmap_meta,

  annotation_colors =
    annotation_colors,

  border_color =
    NA,

  clustering_method =
    "ward.D2",

  filename =
    file.path(
      figdir,
      "Figure2G_sample_transcriptome_dendrogram.png"
    ),

  width = 12,
  height = 8
)


pheatmap(

  heatmap_matrix,

  cluster_rows =
    TRUE,

  cluster_cols =
    TRUE,

  show_rownames =
    FALSE,

  show_colnames =
    FALSE,

  annotation_col =
    heatmap_meta,

  annotation_colors =
    annotation_colors,

  border_color =
    NA,

  clustering_method =
    "ward.D2",

  filename =
    file.path(
      figdir,
      "Figure2G_sample_transcriptome_dendrogram.pdf"
    ),

  width = 12,
  height = 8
)


# ============================================================
# SAVE INDIVIDUAL FIGURE 2 PANELS
# ============================================================

plots_figure2 <- list(

  Figure2A_recurrence =
    p2A,

  Figure2B_disease_breadth =
    p2B,

  Figure2C_isoforms_per_gene =
    p2C,

  Figure2D_length =
    p2D,

  Figure2E_exons =
    p2E,

  Figure2F_category_by_disease =
    p2F
)


for (
  nm in names(
    plots_figure2
  )
) {

  cat(
    "Saving ",
    nm,
    "\n",
    sep = ""
  )


  ggsave(

    filename =
      file.path(
        figdir,
        paste0(
          nm,
          ".png"
        )
      ),

    plot =
      plots_figure2[[nm]],

    width =
      ifelse(
        nm %in%
          c(
            "Figure2D_length",
            "Figure2E_exons"
          ),
        15,
        7
      ),

    height =
      ifelse(
        nm ==
          "Figure2D_length",
        9,
        5.5
      ),

    dpi = 600
  )


  ggsave(

    filename =
      file.path(
        figdir,
        paste0(
          nm,
          ".pdf"
        )
      ),

    plot =
      plots_figure2[[nm]],

    width =
      ifelse(
        nm %in%
          c(
            "Figure2D_length",
            "Figure2E_exons"
          ),
        15,
        7
      ),

    height =
      ifelse(
        nm ==
          "Figure2D_length",
        9,
        5.5
      )
  )
}


# ============================================================
# ASSEMBLE MAIN FIGURE 2
#
# D and E are internally multi-panel, so keep them wide.
# ============================================================

figure2_main <- (

  (
    p2A |
      p2B |
      p2C
  ) /

    p2D /

    p2E /

    p2F

) +

  plot_layout(
    heights =
      c(
        1,
        1.8,
        1,
        1
      )
  )


ggsave(

  filename =
    file.path(
      figdir,
      "Figure2_CWOW_full.png"
    ),

  plot =
    figure2_main,

  width = 18,
  height = 22,
  dpi = 600
)


ggsave(

  filename =
    file.path(
      figdir,
      "Figure2_CWOW_full.pdf"
    ),

  plot =
    figure2_main,

  width = 18,
  height = 22
)


# ============================================================
# SAVE FIGURE 1 C/D/E COMPOSITE
#
# A/B UpSets are saved separately because UpSetR uses grid
# graphics rather than ordinary ggplot objects.
# ============================================================

figure1_CDE <- (

  p1C |
    p1D

) /

  p1E +

  plot_layout(
    heights =
      c(
        1,
        1.4
      )
  )


ggsave(

  filename =
    file.path(
      figdir,
      "Figure1_CDE.png"
    ),

  plot =
    figure1_CDE,

  width = 15,
  height = 12,
  dpi = 600
)


ggsave(

  filename =
    file.path(
      figdir,
      "Figure1_CDE.pdf"
    ),

  plot =
    figure1_CDE,

  width = 15,
  height = 12
)


# ============================================================
# SAVE SOURCE DATA
# ============================================================

fwrite(
  tx_combo,
  file.path(
    figdir,
    "Figure1A_transcript_intersections.tsv"
  ),
  sep = "\t"
)


fwrite(
  gene_combo,
  file.path(
    figdir,
    "Figure1B_gene_intersections.tsv"
  ),
  sep = "\t"
)


fwrite(
  sample_tx_counts,
  file.path(
    figdir,
    "Figure1C_transcript_counts_by_sample.tsv"
  ),
  sep = "\t"
)


fwrite(
  category_summary,
  file.path(
    figdir,
    "Figure1D_category_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  rarefaction_summary,
  file.path(
    figdir,
    "Figure1E_accumulation_curves.tsv"
  ),
  sep = "\t"
)


fwrite(
  support_summary,
  file.path(
    figdir,
    "Figure2A_transcript_recurrence.tsv"
  ),
  sep = "\t"
)


fwrite(
  group_support_summary,
  file.path(
    figdir,
    "Figure2B_disease_group_breadth.tsv"
  ),
  sep = "\t"
)


fwrite(
  gene_tx_by_group,
  file.path(
    figdir,
    "Figure2C_isoforms_per_gene.tsv.gz"
  ),
  sep = "\t"
)


fwrite(
  exon_bar,
  file.path(
    figdir,
    "Figure2E_exon_bins.tsv"
  ),
  sep = "\t"
)


fwrite(
  category_by_group,
  file.path(
    figdir,
    "Figure2F_category_by_disease.tsv"
  ),
  sep = "\t"
)


# ============================================================
# ADDITIONAL CATEGORY SUMMARY BY DISEASE
# ============================================================

category_group_counts <- tx_group_annotation[
  ,
  .(
    N_transcripts =
      uniqueN(
        transcript_id
      )
  ),
  by = .(
    disease_plot,
    category
  )
]


category_group_counts[
  ,
  percent :=
    100 *
    N_transcripts /
    sum(
      N_transcripts
    ),
  by = disease_plot
]


fwrite(
  category_group_counts,
  file.path(
    figdir,
    "CWOW_structural_category_by_disease.tsv"
  ),
  sep = "\t"
)


# ============================================================
# FINAL CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("=============================================\n")
cat("FIGURES COMPLETE\n")
cat("=============================================\n\n")


cat(
  "Output directory:\n",
  figdir,
  "\n\n"
)


cat(
  "Figure 1:\n",
  "  A = transcript UpSet\n",
  "  B = gene UpSet\n",
  "  C = transcript models / individual + P-values\n",
  "  D = structural category composition\n",
  "  E = disease-specific accumulation curves\n\n",
  sep = ""
)


cat(
  "Figure 2:\n",
  "  A = recurrence by structural category\n",
  "  B = disease-group breadth by category\n",
  "  C = transcript models / gene + P-values\n",
  "  D = transcript lengths\n",
  "  E = exon counts\n",
  "  F = category composition by disease\n",
  "  G = transcriptome heatmap/dendrogram\n\n",
  sep = ""
)


cat(
  "Disease colors:\n"
)

print(
  disease_colors
)


cat(
  "\nStructural category colors:\n"
)

print(
  category_colors
)


cat(
  "\nDONE.\n"
)