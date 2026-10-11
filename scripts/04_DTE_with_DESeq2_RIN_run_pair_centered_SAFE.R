# CWOW long-read DTE: DESeq2 
# Source this script in RStudio. No command-line arguments required.
# Input: IsoMatch aggregated full-length (FL) support counts (integers), NOT Oarfish estimates.
# Run all 10 pairwise disease-group comparisons and optional TDP-43 status analysis.
# Design: ~ 0 + group + Age_centered + RIN + run_pair + sex.
# Age_centered is in decades; RIN is centered within each comparison.
# Original workbook RIN and run-pair labels remain available for auditing.
# RIN is exclusively from the original RNA workbook, via augmented QC metadata.
# run_pair is the same sample-level paired-run descriptor used in the batch QC.
# Pairwise filter: >=5 FL in >=5 samples across the two comparison groups.
# New outputs are written to results/DTE_DESeq2/moderate/.

# ---------------------------- USER SETTINGS -----------------------------
BASE <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs"
PROJECT <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts"
SCENARIO <- "moderate"
COUNT_FILE <- file.path(BASE, "04_confidence_filter_rescuedLR_274", SCENARIO,
                        "CWOW_274_rescuedLR.confidence_filtered.counts.tsv.gz")
MAP_FILE <- file.path(BASE, "04_confidence_filter_rescuedLR_274", SCENARIO,
                      "CWOW_274_rescuedLR.transcript_to_gene.tsv.gz")
METADATA_FILE <- file.path(dirname(BASE), "QC_technical_RIN_normalization_CCA",
                           "CWOW_n274_DTE_metadata_RIN_SMRTcell_RGID_augmented.tsv")
OUT <- file.path(dirname(BASE), "DTE_DESeq2", SCENARIO)

FDR <- 0.05
MIN_FL_PER_SAMPLE <- 5L
MIN_SAMPLES_WITH_FL <- 5L
MIN_GROUP_N <- 3L
RUN_COMPARISONS <- "ALL" # all 10 disease-group comparisons
SKIP_EXISTING_DTE <- FALSE # new design: do not reuse older-model DTE outputs
# ------------------------------------------------------------------------

if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table")
if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install DESeq2 with BiocManager")
library(data.table)
for (p in c(COUNT_FILE, MAP_FILE, METADATA_FILE)) {
  if (!file.exists(p)) stop("Missing input: ", p)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

comparisons <- list(
  AD_vs_CONTROL = c("AD", "CONTROL"),
  LBD_S_vs_CONTROL = c("LBD_S", "CONTROL"),
  LBD_AS_vs_CONTROL = c("LBD_AS", "CONTROL"),
  LBD_ATS_vs_CONTROL = c("LBD_ATS", "CONTROL"),
  LBD_S_vs_AD = c("LBD_S", "AD"),
  LBD_AS_vs_AD = c("LBD_AS", "AD"),
  LBD_ATS_vs_AD = c("LBD_ATS", "AD"),
  LBD_AS_vs_LBD_S = c("LBD_AS", "LBD_S"),
  LBD_ATS_vs_LBD_S = c("LBD_ATS", "LBD_S"),
  LBD_ATS_vs_LBD_AS = c("LBD_ATS", "LBD_AS")
)
RUN_TDP43 <- TRUE  # Additional TDP.43-positive vs negative comparison

if (identical(RUN_COMPARISONS, "ALL")) {
  selected <- names(comparisons)
} else {
  selected <- RUN_COMPARISONS
  if (!length(selected) || !all(selected %in% names(comparisons)))
    stop("RUN_COMPARISONS must be 'ALL' or one/more names from comparisons")
}

message("Reading corrected QC metadata (workbook RIN and JSON-derived runs)...")
meta <- fread(METADATA_FILE, check.names = FALSE)
required <- c("sample_id", "group", "sex", "Age", "RIN", "RIN_workbook",
              "all_pacbio_runs", "rgid_1", "rgid_2", "TDP.43")
if (!all(required %in% names(meta)))
  stop("Missing augmented metadata columns: ", paste(setdiff(required, names(meta)), collapse = ", "))
meta[, sample_id := trimws(as.character(sample_id))]
meta[, group := gsub("-", "_", toupper(trimws(as.character(group))))]
meta[, sex := toupper(trimws(as.character(sex)))]
meta[sex == "MALE", sex := "XY"]
meta[sex == "FEMALE", sex := "XX"]
meta[, Age := suppressWarnings(as.numeric(Age))]
meta[, RIN := suppressWarnings(as.numeric(RIN_workbook))]
meta[, TDP43 := suppressWarnings(as.numeric(`TDP.43`))]
if (any(!is.na(meta$TDP43) & !meta$TDP43 %in% c(0, 1)))
  stop("TDP.43 must be 0, 1 or missing")
# Same run-pair definition as 04_CWOW_SMRT_Run_Confounding_PC_Assessment.R.
# Each sample contributes one combined count profile, not two independent replicates.
clean_pair <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(NA_character_)
  tokens <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  tokens <- tokens[nzchar(tokens)]
  if (!length(tokens)) return(NA_character_)
  paste(sort(unique(tokens)), collapse = ";")
}
meta[, run_pair_original := vapply(as.character(all_pacbio_runs), clean_pair, character(1))]
# Stable, R-safe factor labels. Preserve the exact original run pair for audit.
run_keys <- sort(unique(meta$run_pair_original[!is.na(meta$run_pair_original)]))
run_map <- data.table(run_pair_original = run_keys,
                      run_pair = sprintf("RunPair_%02d", seq_along(run_keys)))
meta[run_map, on = "run_pair_original", run_pair := i.run_pair]

if (anyDuplicated(meta$sample_id)) stop("Duplicate metadata sample IDs")
if (any(!meta$sex %in% c("XX", "XY"))) stop("Unexpected sex labels")
if (any(!meta$group %in% c("CONTROL", "AD", "LBD_S", "LBD_AS", "LBD_ATS")))
  stop("Unexpected disease group labels")
if (any(is.na(meta$run_pair))) stop("Missing JSON-derived run pair: fix metadata before DTE")
if (any(is.na(meta$rgid_1) | is.na(meta$rgid_2))) stop("Missing SMRT read groups")

fwrite(run_map, file.path(OUT, "run_pair_safe_label_mapping.tsv"), sep = "\t")
message("Reading count matrix (large file)...")
raw <- fread(COUNT_FILE, check.names = FALSE, showProgress = TRUE)
if (!identical(names(raw)[1], "transcript_id")) stop("First column must be transcript_id")
ids <- raw$transcript_id
if (anyDuplicated(ids)) stop("Duplicate transcript IDs")
samples <- names(raw)[-1L]
if (length(samples) != 274L) stop("Expected 274 count samples; found ", length(samples))
if (anyDuplicated(samples)) stop("Duplicate count sample names")
missing_meta <- setdiff(samples, meta$sample_id)
if (length(missing_meta)) stop("Missing metadata for: ", paste(head(missing_meta, 20), collapse = ", "))
meta <- meta[match(samples, sample_id)]
stopifnot(identical(meta$sample_id, samples))
if (anyNA(meta[, .(Age, sex, group, run_pair)])) stop("Missing required model covariates")
# The workbook supplies 273/274 RIN values; omit the missing-RIN individual
# from each comparison, never impute or substitute Tapestation_DIN_or_RIN.
message("Samples with workbook RIN: ", sum(!is.na(meta$RIN)), " / ", nrow(meta))
fwrite(meta[is.na(RIN), .(sample_id, group, RIN_workbook)],
       file.path(OUT, "samples_excluded_missing_workbook_RIN.tsv"), sep = "\t")

message("Reading transcript-to-gene mapping...")
map <- fread(MAP_FILE, select = c("transcript_id", "gene_id"))
if (anyDuplicated(map$transcript_id)) stop("Duplicate transcript-to-gene mappings")
map <- map[match(ids, transcript_id)]
if (anyNA(map$gene_id)) stop("Some transcripts lack gene assignments")

# Build one dense integer count matrix; filter per comparison before DESeqDataSet construction.
# Do not apply a weak whole-cohort filter first, and avoid retaining extra full-sized matrices.
message("Converting counts to matrix...")
count_cols <- names(raw)[-1L]
counts <- as.matrix(raw[, ..count_cols])
rownames(counts) <- ids
rm(raw)
invisible(gc())
if (anyNA(counts) || any(counts < 0) || any(counts != round(counts)))
  stop("Counts must be nonnegative integers without missing values")
if (any(counts > .Machine$integer.max)) stop("Count exceeds R integer range")
storage.mode(counts) <- "integer"

meta_export <- copy(meta)
meta_export[, Age_centered := (Age - mean(Age, na.rm = TRUE)) / 10]
# The RIN column here remains the raw workbook RIN; per-comparison centered
# values actually used in fitting are saved to model_metadata_used.tsv.
fwrite(meta_export, file.path(OUT, "matched_metadata_all_original_columns.tsv"), sep = "\t")
fwrite(as.data.table(table(meta$group, meta$sex)), file.path(OUT, "group_by_sex.tsv"), sep = "\t")

# Exact user-requested model; group coefficients are group-specific adjusted means.
# Validate full rank independently for each pairwise contrast after removing
# samples with missing workbook RIN and unused factor levels.
make_design <- function(md, tdp43 = FALSE) {
  f <- if (tdp43) {
    ~ 0 + group + Age_centered + RIN + run_pair + sex + TDP43_status
  } else {
    ~ 0 + group + Age_centered + RIN + run_pair + sex
  }
  mm <- model.matrix(f, data = md)
  if (qr(mm)$rank < ncol(mm)) {
    stop("Rank-deficient design: ", deparse(f),
         "; check run-pair / disease / sex confounding for this subset")
  }
  if (nrow(mm) <= ncol(mm)) stop("Insufficient residual degrees of freedom")
  f
}

library(DESeq2)
filter_rows <- list()
record_filter <- function(nm, group1, group2, n1, n2, n_input, n_pass, analysis) {
  filter_rows[[nm]] <<- data.table(
    analysis = analysis, comparison = nm,
    group_1 = group1, n_group_1 = as.integer(n1),
    group_2 = group2, n_group_2 = as.integer(n2),
    n_total = as.integer(n1 + n2),
    n_transcripts_input = as.integer(n_input),
    n_transcripts_passing = as.integer(n_pass),
    percent_retained = round(100 * n_pass / n_input, 2),
    min_fl_per_sample = MIN_FL_PER_SAMPLE,
    min_samples_with_fl = MIN_SAMPLES_WITH_FL)
  # Written together at the end, in a stable comparison order.
}
for (nm in selected) {
  message("\n========== DTE: ", nm, " ==========")
  pair <- comparisons[[nm]]
  ix <- which(meta$group %in% pair & !is.na(meta$RIN))
  md <- as.data.frame(meta[ix])
  rownames(md) <- md$sample_id
  md$group <- factor(md$group, levels = c(pair[2], pair[1]))
  md$run_pair <- droplevels(factor(md$run_pair))
  md$sex <- droplevels(factor(md$sex))
  # Keep the requested formula, but center/scale the numeric covariates.
  # RIN values are sourced ONLY from RIN_workbook above (never imputed).
  # Centering RIN changes the coefficient reference point, not group contrasts.
  md$Age_centered <- (md$Age - mean(md$Age)) / 10  # per 10 years
  md$RIN <- md$RIN - mean(md$RIN)                 # centered RIN units
  if (sd(md$Age_centered) == 0 || sd(md$RIN) == 0)
    stop("Age or workbook RIN has no variation in ", nm)
  group_n <- table(md$group)
  if (length(group_n) != 2L || any(group_n < MIN_GROUP_N))
    stop("Insufficient samples in ", nm, ": ", paste(names(group_n), group_n, collapse = ", "))
  f <- make_design(md)
  out <- file.path(OUT, nm)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(deparse(f),
               "Age_centered = (Age - mean(Age within comparison))/10",
               "RIN = RIN_workbook - mean(RIN_workbook within comparison)",
               "run_pair uses safe IDs; see run_pair_safe_label_mapping.tsv"),
             file.path(out, "model_formula.txt"))
  fwrite(as.data.table(md), file.path(out, "model_metadata_used.tsv"), sep = "\t")
  fwrite(data.table(group = names(group_n), n = as.integer(group_n)),
         file.path(out, "group_sample_sizes.tsv"), sep = "\t")

  # Prevalence is assessed across BOTH groups combined: a transcript can qualify
  # even if all qualifying samples are from one disease group.
  # A threshold of >=5 FL reads in >=5 samples is an exploratory choice.
  sub <- counts[, ix, drop = FALSE]
  keep <- rowSums(sub >= MIN_FL_PER_SAMPLE) >= MIN_SAMPLES_WITH_FL
  message("  Retained ", sum(keep), " / ", length(keep), " transcripts with >=",
          MIN_FL_PER_SAMPLE, " FL in >=", MIN_SAMPLES_WITH_FL, " samples")
  fwrite(data.table(comparison = nm, n_samples = ncol(sub), n_transcripts_input = nrow(sub),
                    n_transcripts_tested = sum(keep), min_fl = MIN_FL_PER_SAMPLE,
                    min_samples = MIN_SAMPLES_WITH_FL),
         file.path(out, "filter_summary.tsv"), sep = "\t")
  record_filter(nm, pair[1], pair[2], unname(group_n[pair[1]]),
                unname(group_n[pair[2]]), nrow(sub), sum(keep), "Disease_pairwise")
  if (!any(keep)) stop("No transcripts pass expression filter for ", nm)
  if (SKIP_EXISTING_DTE && all(file.exists(file.path(out, c(
      "DTE_results.tsv.gz", "DTE_FDR05.tsv.gz", "size_factors.tsv", "run_summary.rds"))))) {
    message("  Existing DTE outputs found; skipping DESeq2 fit (filtering statistics refreshed).")
    rm(sub, keep)
    invisible(gc())
    next
  }

  dds <- DESeqDataSetFromMatrix(countData = sub[keep, , drop = FALSE],
                               colData = md, design = f)
  rm(sub)
  invisible(gc())
  dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)
  # With ~0 + group, use a numeric contrast: group1 - group2.
  coef_names <- resultsNames(dds)
  g1 <- paste0("group", pair[1]); g2 <- paste0("group", pair[2])
  if (!all(c(g1, g2) %in% coef_names))
    stop("Group coefficients absent: ", paste(setdiff(c(g1, g2), coef_names), collapse = ", "))
  contrast_vector <- setNames(numeric(length(coef_names)), coef_names)
  contrast_vector[g1] <- 1; contrast_vector[g2] <- -1
  res <- results(dds, contrast = unname(contrast_vector), alpha = FDR)
  ans <- as.data.table(as.data.frame(res), keep.rownames = "transcript_id")
  ans[, gene_id := map$gene_id[match(transcript_id, ids)]]
  setcolorder(ans, c("transcript_id", "gene_id", setdiff(names(ans), c("transcript_id", "gene_id"))))
  setorder(ans, padj)
  fwrite(ans, file.path(out, "DTE_results.tsv.gz"), sep = "\t")
  fwrite(ans[!is.na(padj) & padj < FDR], file.path(out, "DTE_FDR05.tsv.gz"), sep = "\t")
  fwrite(data.table(sample_id = colnames(dds), size_factor = sizeFactors(dds)),
         file.path(out, "size_factors.tsv"), sep = "\t")
  saveRDS(list(comparison = nm, design = deparse(f), n_tested = nrow(ans),
               n_fdr05 = sum(!is.na(ans$padj) & ans$padj < FDR),
               session = sessionInfo()), file.path(out, "run_summary.rds"))
  fwrite(data.table(comparison = nm, n_initial_pass = nrow(ans),
                    n_padj_nonNA = sum(!is.na(ans$padj)),
                    n_padj_NA = sum(is.na(ans$padj)),
                    n_significant = sum(!is.na(ans$padj) & ans$padj < FDR)),
         file.path(out, "DESeq2_independent_filtering_summary.tsv"), sep = "\t")
  message("  Significant transcripts (FDR < ", FDR, "): ",
          sum(!is.na(ans$padj) & ans$padj < FDR))
  rm(dds, res, ans, keep)
  invisible(gc())
}
# TDP-43 comparison uses only classified samples; disease group is a covariate.
# Missing TDP.43 classifications are NOT interpreted as negative.
if (RUN_TDP43) {
  nm <- "TDP43_positive_vs_negative"
  message("\n========== DTE: ", nm, " ==========")
  ix <- which(!is.na(meta$TDP43) & !is.na(meta$RIN))
  md <- as.data.frame(meta[ix])
  rownames(md) <- md$sample_id
  md$TDP43_status <- factor(ifelse(md$TDP43 == 1, "Positive", "Negative"),
                            levels = c("Negative", "Positive"))
  md$group <- droplevels(factor(md$group))
  md$run_pair <- droplevels(factor(md$run_pair))
  md$sex <- droplevels(factor(md$sex))
  # Keep the requested formula, but center/scale the numeric covariates.
  # RIN values are sourced ONLY from RIN_workbook above (never imputed).
  # Centering RIN changes the coefficient reference point, not group contrasts.
  md$Age_centered <- (md$Age - mean(md$Age)) / 10  # per 10 years
  md$RIN <- md$RIN - mean(md$RIN)                 # centered RIN units
  if (sd(md$Age_centered) == 0 || sd(md$RIN) == 0)
    stop("Age or workbook RIN has no variation in ", nm)
  status_n <- table(md$TDP43_status)
  if (any(status_n < MIN_GROUP_N)) stop("Too few TDP-43 classified samples")
  f <- make_design(md, tdp43 = TRUE)
  out <- file.path(OUT, nm)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(deparse(f),
               "Age_centered = (Age - mean(Age within comparison))/10",
               "RIN = RIN_workbook - mean(RIN_workbook within comparison)",
               "run_pair uses safe IDs; see run_pair_safe_label_mapping.tsv"),
             file.path(out, "model_formula.txt"))
  fwrite(as.data.table(md), file.path(out, "model_metadata_used.tsv"), sep = "\t")
  fwrite(data.table(status = names(status_n), n = as.integer(status_n)),
         file.path(out, "group_sample_sizes.tsv"), sep = "\t")
  fwrite(as.data.table(table(md$group, md$TDP43_status)),
         file.path(out, "neuropathology_by_TDP43.tsv"), sep = "\t")
  sub <- counts[, ix, drop = FALSE]
  keep <- rowSums(sub >= MIN_FL_PER_SAMPLE) >= MIN_SAMPLES_WITH_FL
  message("  Retained ", sum(keep), " / ", length(keep), " transcripts")
  record_filter(nm, "TDP43_Positive", "TDP43_Negative",
                status_n[["Positive"]], status_n[["Negative"]],
                nrow(sub), sum(keep), "TDP43_adjusted_for_neuropathology")
  fwrite(data.table(comparison = nm, n_samples = ncol(sub),
                    n_transcripts_input = nrow(sub), n_transcripts_tested = sum(keep),
                    min_fl = MIN_FL_PER_SAMPLE, min_samples = MIN_SAMPLES_WITH_FL),
         file.path(out, "filter_summary.tsv"), sep = "\t")
  if (!any(keep)) stop("No transcripts pass TDP-43 filter")
  if (SKIP_EXISTING_DTE && all(file.exists(file.path(out, c(
      "DTE_results.tsv.gz", "DTE_FDR05.tsv.gz", "size_factors.tsv", "run_summary.rds"))))) {
    message("  Existing TDP-43 DTE outputs found; skipping DESeq2 fit.")
  } else {
  dds <- DESeqDataSetFromMatrix(countData = sub[keep, , drop = FALSE],
                                 colData = md, design = f)
  rm(sub)
  invisible(gc())
  dds <- DESeq(dds, quiet = TRUE, parallel = FALSE)
  res <- results(dds, contrast = c("TDP43_status", "Positive", "Negative"), alpha = FDR)
  ans <- as.data.table(as.data.frame(res), keep.rownames = "transcript_id")
  ans[, gene_id := map$gene_id[match(transcript_id, ids)]]
  setcolorder(ans, c("transcript_id", "gene_id",
                     setdiff(names(ans), c("transcript_id", "gene_id"))))
  setorder(ans, padj)
  fwrite(ans, file.path(out, "DTE_results.tsv.gz"), sep = "\t")
  fwrite(ans[!is.na(padj) & padj < FDR],
         file.path(out, "DTE_FDR05.tsv.gz"), sep = "\t")
  fwrite(data.table(sample_id = colnames(dds), size_factor = sizeFactors(dds)),
         file.path(out, "size_factors.tsv"), sep = "\t")
  saveRDS(list(comparison = nm, design = deparse(f), n_tested = nrow(ans),
               n_fdr05 = sum(!is.na(ans$padj) & ans$padj < FDR),
               session = sessionInfo()), file.path(out, "run_summary.rds"))
  fwrite(data.table(comparison = nm, n_initial_pass = nrow(ans),
                    n_padj_nonNA = sum(!is.na(ans$padj)),
                    n_padj_NA = sum(is.na(ans$padj)),
                    n_significant = sum(!is.na(ans$padj) & ans$padj < FDR)),
         file.path(out, "DESeq2_independent_filtering_summary.tsv"), sep = "\t")
  message("  Significant transcripts (FDR < ", FDR, "): ",
          sum(!is.na(ans$padj) & ans$padj < FDR))
  rm(dds, res, ans, keep)
  invisible(gc())
  } # end fit-or-skip for TDP-43
  if (exists("sub", inherits = FALSE)) rm(sub)
  if (exists("keep", inherits = FALSE)) rm(keep)
}

# One complete supplementary table: 10 disease comparisons plus TDP-43.
filter_table <- rbindlist(filter_rows, use.names = TRUE, fill = TRUE)
filter_order <- c(names(comparisons), "TDP43_positive_vs_negative")
filter_table[, sort_order := match(comparison, filter_order)]
setorder(filter_table, sort_order)
filter_table[, sort_order := NULL]
fwrite(filter_table, file.path(OUT, "Supplementary_Table_DTE_Transcript_Filtering.tsv"), sep = "\t")
message("Saved filtering table with ", nrow(filter_table), " comparisons.")
message("\nDTE completed. Results: ", OUT)
message("Supplementary table: ",
        file.path(OUT, "Supplementary_Table_DTE_Transcript_Filtering.tsv"))
