# ============================================================
# CWOW PACBIO KINNEX
# HIGH-CONFIDENCE COHORT TRANSCRIPTOME SUMMARY
#
# CURRENT COHORT:
#   211 completed samples
#
# Calculates:
#
#   1. Number of sample-level transcript observations BEFORE
#      IsoMatch cohort merging
#
#   2. Number of UNIQUE high-confidence transcript structures
#      AFTER IsoMatch harmonization
#
#   3. Number of unique genes represented
#
#   4. FSM / NIC / NNC composition of the UNIQUE harmonized
#      high-confidence transcriptome
#
#   5. Number of individuals supporting each cohort transcript
#
#   6. Number of transcripts supported by >= 1, 2, 3, 5, 10,
#      20, and 50 individuals
#
#   7. Abstract-ready summary sentences
#
# ============================================================


# ============================================================
# 0. LIBRARIES
# ============================================================

library(data.table)


# ============================================================
# 1. PROJECT PATHS
# ============================================================

outdir <- paste0(
  "/tgen_labs/jfryer/kolney/LBD_CWOW/",
  "isoform_CWOW/isomatch_cohort"
)


# Final high-confidence transcript table generated previously
#
# One row per UNIQUE high-confidence ISOMT transcript
#
hc_file <- file.path(
  outdir,
  "CWOW_high_confidence_transcripts.tsv"
)


# IsoMatch source-to-merged transcript mapping
#
# Each row represents a sample-level PB transcript contributing
# to a merged ISOMT transcript.
#
track_file <- file.path(
  outdir,
  "CWOW_cohort.track.tsv.gz"
)


# Full IsoMatch classification of the harmonized cohort
#
classification_file <- file.path(
  outdir,
  "CWOW_cohort_gencodev38.classification.txt.gz"
)


# ============================================================
# 2. CHECK INPUT FILES
# ============================================================

required_files <- c(
  hc_file,
  track_file,
  classification_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat(
    "\nMissing required files:\n"
  )

  print(
    missing_files
  )

  stop(
    "Required input files are missing."
  )
}


cat(
  "\n=============================================\n"
)

cat(
  "CWOW HIGH-CONFIDENCE TRANSCRIPTOME SUMMARY\n"
)

cat(
  "=============================================\n\n"
)


# ============================================================
# 3. READ HIGH-CONFIDENCE TRANSCRIPT TABLE
# ============================================================

cat(
  "Reading high-confidence cohort transcript table...\n"
)

hc <- fread(
  hc_file
)

cat(
  "Rows in HC table:",
  format(
    nrow(hc),
    big.mark = ","
  ),
  "\n"
)


# Make sure transcript IDs are unique
stopifnot(
  uniqueN(
    hc$transcript_id
  ) ==
    nrow(hc)
)


# ============================================================
# 4. NUMBER OF UNIQUE COHORT TRANSCRIPT STRUCTURES
#
# THIS is the harmonized cohort-level transcript count.
# ============================================================

N_unique_transcripts <- uniqueN(
  hc$transcript_id
)


cat(
  "\nUnique high-confidence cohort transcripts:",
  format(
    N_unique_transcripts,
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# 5. READ ISOMATCH TRACK FILE
#
# This maps:
#
#   sample PB transcript
#
#           ->
#
#   cohort ISOMT transcript
#
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
    "merged_gene_id",
    "src_tx_id",
    "src_gene_id",
    "src_file_name"
  )
)


setnames(
  track,
  old = c(
    "merged_tx_id",
    "merged_gene_id",
    "src_tx_id",
    "src_gene_id",
    "src_file_name"
  ),
  new = c(
    "transcript_id",
    "isom_gene_id",
    "pbid",
    "pb_gene_id",
    "src_file_name"
  )
)


# Remove filename suffix to create sample name
track[
  ,
  sample :=
    sub(
      "\\.filtered\\.gtf$",
      "",
      src_file_name
    )
]


cat(
  "Samples represented:",
  uniqueN(
    track$sample
  ),
  "\n"
)


# ============================================================
# 6. RESTRICT TRACK TO FINAL HC TRANSCRIPTS
#
# This is critical.
#
# We want the number of sample-level transcript observations
# that contributed specifically to the FINAL 1.807M HC catalog.
# ============================================================

track_hc <- track[
  transcript_id %chin%
    hc$transcript_id
]


cat(
  "\nSource PB transcript mappings contributing to HC catalog:",
  format(
    nrow(track_hc),
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# 7. NUMBER OF TRANSCRIPT OBSERVATIONS BEFORE COHORT MERGING
#
# Each unique:
#
#   sample + PB transcript
#
# represents one sample-level transcript observation.
#
# We explicitly count sample + pbid because PB IDs can repeat
# between samples.
# ============================================================

N_transcript_observations <- uniqueN(
  track_hc[
    ,
    paste(
      sample,
      pbid,
      sep = "::"
    )
  ]
)


cat(
  "\nSample-level HC transcript observations BEFORE merging:",
  format(
    N_transcript_observations,
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# 8. CHECK FOR DUPLICATE SAMPLE/PB PAIRS
# ============================================================

duplicate_source_ids <- track_hc[
  ,
  .N,
  by = .(
    sample,
    pbid
  )
][
  N > 1
]


cat(
  "Duplicated sample/PB transcript combinations:",
  nrow(
    duplicate_source_ids
  ),
  "\n"
)


# ============================================================
# 9. COLLAPSE RATIO
#
# How many sample-level transcript observations correspond
# to each unique cohort model, on average?
# ============================================================

observation_to_unique_ratio <-
  N_transcript_observations /
  N_unique_transcripts


cat(
  "\nObservation-to-unique-transcript ratio:",
  round(
    observation_to_unique_ratio,
    2
  ),
  "\n"
)


# ============================================================
# 10. UNIQUE ISOMATCH GENES / LOCI
#
# This is the number of cohort gene clusters produced by
# IsoMatch.
#
# This includes novel loci and therefore is different from
# counting only GENCODE gene IDs.
# ============================================================

hc_gene_map <- unique(
  track_hc[
    ,
    .(
      transcript_id,
      isom_gene_id
    )
  ]
)


N_isomatch_genes <- uniqueN(
  hc_gene_map$isom_gene_id
)


cat(
  "\nUnique IsoMatch gene loci:",
  format(
    N_isomatch_genes,
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# 11. READ ISOMATCH CLASSIFICATION
#
# Needed for reference gene assignments and structural classes.
# ============================================================

cat(
  "\nReading cohort classification...\n"
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
    "structural_category",
    "subcategory",
    "ref_gene_id",
    "ref_gene_name",
    "ref_tx_id"
  )
)


setnames(
  classification,
  "isoform_id",
  "transcript_id"
)


# Keep only final HC transcripts
classification_hc <- classification[
  transcript_id %chin%
    hc$transcript_id
]


# ============================================================
# 12. VERIFY HC CLASSIFICATION SIZE
# ============================================================

cat(
  "\nHC transcripts with IsoMatch classification:",
  format(
    nrow(
      classification_hc
    ),
    big.mark = ","
  ),
  "\n"
)


stopifnot(
  nrow(
    classification_hc
  ) ==
    N_unique_transcripts
)


# ============================================================
# 13. UNIQUE GENCODE GENES REPRESENTED
#
# This counts only transcripts with a GENCODE gene assignment.
#
# Novel/intergenic transcripts without reference gene
# assignments are not counted here.
# ============================================================

N_reference_genes <- classification_hc[
  !is.na(
    ref_gene_id
  ) &
    ref_gene_id != "",
  uniqueN(
    ref_gene_id
  )
]


cat(
  "Unique GENCODE genes represented:",
  format(
    N_reference_genes,
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# 14. STRUCTURAL CATEGORY SUMMARY
#
# FSM / NIC / NNC
# ============================================================

classification_summary <- classification_hc[
  ,
  .(
    N = .N
  ),
  by = structural_category
]


classification_summary[
  ,
  percent :=
    100 *
      N /
      sum(N)
]


setorder(
  classification_summary,
  -N
)


cat(
  "\n=============================================\n"
)

cat(
  "STRUCTURAL CATEGORY SUMMARY\n"
)

cat(
  "=============================================\n\n"
)


print(
  classification_summary
)


# ============================================================
# 15. EXTRACT FSM / NIC / NNC VALUES
# ============================================================

N_FSM <- classification_summary[
  structural_category ==
    "full-splice_match",
  N
]


N_NIC <- classification_summary[
  structural_category ==
    "novel_in_catalog",
  N
]


N_NNC <- classification_summary[
  structural_category ==
    "novel_not_in_catalog",
  N
]


PCT_FSM <- classification_summary[
  structural_category ==
    "full-splice_match",
  percent
]


PCT_NIC <- classification_summary[
  structural_category ==
    "novel_in_catalog",
  percent
]


PCT_NNC <- classification_summary[
  structural_category ==
    "novel_not_in_catalog",
  percent
]


# ============================================================
# 16. VERIFY THAT HC CATALOG CONTAINS ONLY FSM / NIC / NNC
# ============================================================

expected_categories <- c(
  "full-splice_match",
  "novel_in_catalog",
  "novel_not_in_catalog"
)


unexpected_categories <- classification_summary[
  !structural_category %chin%
    expected_categories
]


if (
  nrow(
    unexpected_categories
  ) > 0
) {

  warning(
    "Unexpected structural categories found in HC catalog."
  )

  print(
    unexpected_categories
  )
}


# ============================================================
# 17. SAMPLE SUPPORT FOR EACH UNIQUE COHORT TRANSCRIPT
#
# We calculate this directly from the IsoMatch track file.
# ============================================================

transcript_support <- track_hc[
  ,
  .(
    n_samples =
      uniqueN(
        sample
      ),

    n_source_transcripts =
      uniqueN(
        paste(
          sample,
          pbid,
          sep = "::"
        )
      )
  ),
  by = .(
    transcript_id,
    isom_gene_id
  )
]


cat(
  "\nTranscripts with support information:",
  format(
    nrow(
      transcript_support
    ),
    big.mark = ","
  ),
  "\n"
)


# ============================================================
# 18. SUPPORT THRESHOLD SUMMARY
# ============================================================

support_thresholds <- c(
  1,
  2,
  3,
  5,
  10,
  20,
  50
)


support_summary <- rbindlist(
  lapply(
    support_thresholds,
    function(k) {

      data.table(

        min_individuals =
          k,

        N_transcripts =
          sum(
            transcript_support$n_samples >= k
          )
      )
    }
  )
)


support_summary[
  ,
  percent_of_HC_catalog :=
    100 *
      N_transcripts /
      N_unique_transcripts
]


cat(
  "\n=============================================\n"
)

cat(
  "COHORT SUPPORT SUMMARY\n"
)

cat(
  "=============================================\n\n"
)


print(
  support_summary
)


# ============================================================
# 19. SUPPORT BY STRUCTURAL CATEGORY
# ============================================================

support_with_class <- merge(
  transcript_support,
  classification_hc[
    ,
    .(
      transcript_id,
      structural_category,
      subcategory
    )
  ],
  by = "transcript_id",
  all.x = TRUE
)


support_by_category <- support_with_class[
  ,
  .(
    N_transcripts =
      .N,

    median_samples =
      as.numeric(
        median(
          n_samples
        )
      ),

    mean_samples =
      as.numeric(
        mean(
          n_samples
        )
      ),

    singleton =
      as.numeric(
        sum(
          n_samples == 1
        )
      ),

    pct_singleton =
      as.numeric(
        100 *
          mean(
            n_samples == 1
          )
      ),

    N_2plus =
      as.numeric(
        sum(
          n_samples >= 2
        )
      ),

    N_3plus =
      as.numeric(
        sum(
          n_samples >= 3
        )
      ),

    N_5plus =
      as.numeric(
        sum(
          n_samples >= 5
        )
      ),

    N_10plus =
      as.numeric(
        sum(
          n_samples >= 10
        )
      )
  ),
  by = structural_category
]


cat(
  "\n=============================================\n"
)

cat(
  "SUPPORT BY STRUCTURAL CATEGORY\n"
)

cat(
  "=============================================\n\n"
)


print(
  support_by_category
)


# ============================================================
# 20. NUMBER OF SAMPLES IN CURRENT COHORT
# ============================================================

N_samples <- uniqueN(
  track_hc$sample
)


# ============================================================
# 21. CREATE MASTER SUMMARY TABLE
# ============================================================

master_summary <- data.table(

  metric = c(

    "Samples currently included",

    "Sample-level transcript observations before IsoMatch merging",

    "Unique high-confidence transcript models after IsoMatch harmonization",

    "Unique IsoMatch gene loci",

    "Unique GENCODE genes represented",

    "FSM transcripts",

    "NIC transcripts",

    "NNC transcripts",

    "FSM percent",

    "NIC percent",

    "NNC percent",

    "Mean source observations per unique transcript"
  ),

  value = c(

    N_samples,

    N_transcript_observations,

    N_unique_transcripts,

    N_isomatch_genes,

    N_reference_genes,

    N_FSM,

    N_NIC,

    N_NNC,

    PCT_FSM,

    PCT_NIC,

    PCT_NNC,

    observation_to_unique_ratio
  )
)


cat(
  "\n=============================================\n"
)

cat(
  "MASTER SUMMARY\n"
)

cat(
  "=============================================\n\n"
)


print(
  master_summary
)


# ============================================================
# 22. SAVE OUTPUT TABLES
# ============================================================

fwrite(
  master_summary,
  file.path(
    outdir,
    "CWOW_high_confidence_cohort_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  classification_summary,
  file.path(
    outdir,
    "CWOW_HC_FSM_NIC_NNC_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  transcript_support,
  file.path(
    outdir,
    "CWOW_HC_transcript_sample_support.tsv.gz"
  ),
  sep = "\t"
)


fwrite(
  support_summary,
  file.path(
    outdir,
    "CWOW_HC_support_threshold_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  support_by_category,
  file.path(
    outdir,
    "CWOW_HC_support_by_category.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 23. PRINT ABSTRACT-READY RESULTS
# ============================================================

cat(
  "\n\n=============================================\n"
)

cat(
  "ABSTRACT-READY NUMBERS\n"
)

cat(
  "=============================================\n\n"
)


cat(
  "Current cohort: ",
  N_samples,
  " individuals\n\n",
  sep = ""
)


cat(
  "Sample-level transcript observations before cohort merging: ",
  format(
    N_transcript_observations,
    big.mark = ","
  ),
  "\n",
  sep = ""
)


cat(
  "Unique high-confidence transcript models after harmonization: ",
  format(
    N_unique_transcripts,
    big.mark = ","
  ),
  "\n",
  sep = ""
)


cat(
  "Unique IsoMatch gene loci: ",
  format(
    N_isomatch_genes,
    big.mark = ","
  ),
  "\n",
  sep = ""
)


cat(
  "Unique GENCODE genes represented: ",
  format(
    N_reference_genes,
    big.mark = ","
  ),
  "\n\n",
  sep = ""
)


cat(
  "FSM: ",
  format(
    N_FSM,
    big.mark = ","
  ),
  " (",
  round(
    PCT_FSM,
    1
  ),
  "%)\n",
  sep = ""
)


cat(
  "NIC: ",
  format(
    N_NIC,
    big.mark = ","
  ),
  " (",
  round(
    PCT_NIC,
    1
  ),
  "%)\n",
  sep = ""
)


cat(
  "NNC: ",
  format(
    N_NNC,
    big.mark = ","
  ),
  " (",
  round(
    PCT_NNC,
    1
  ),
  "%)\n\n",
  sep = ""
)


# ============================================================
# 24. PRINT PROVISIONAL ABSTRACT SENTENCE
# ============================================================

cat(
  "PROVISIONAL ABSTRACT SENTENCE:\n\n"
)


cat(
  paste0(
    "Across ",
    N_samples,
    " individuals, ",
    format(
      N_transcript_observations,
      big.mark = ","
    ),
    " high-confidence sample-level transcript observations ",
    "were harmonized into ",
    format(
      N_unique_transcripts,
      big.mark = ","
    ),
    " unique cohort transcript models representing ",
    format(
      N_reference_genes,
      big.mark = ","
    ),
    " annotated genes. The harmonized transcriptome comprised ",
    format(
      N_FSM,
      big.mark = ","
    ),
    " (",
    round(
      PCT_FSM,
      1
    ),
    "%) full-splice matches, ",
    format(
      N_NIC,
      big.mark = ","
    ),
    " (",
    round(
      PCT_NIC,
      1
    ),
    "%) novel-in-catalog transcripts, and ",
    format(
      N_NNC,
      big.mark = ","
    ),
    " (",
    round(
      PCT_NNC,
      1
    ),
    "%) novel-not-in-catalog transcripts.\n"
  )
)


cat(
  "\n=============================================\n"
)

cat(
  "DONE\n"
)

cat(
  "=============================================\n"
)