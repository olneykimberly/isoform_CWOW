# CWOW long-read DTE: DESeq2 
# Source this script in RStudio. No command-line arguments required.
# Input: IsoMatch aggregated full-length (FL) support counts (integers), NOT Oarfish estimates.
# Run all 10 pairwise disease-group comparisons and optional TDP-43 status analysis.

# ---------------------------- USER SETTINGS -----------------------------
BASE <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs"
PROJECT <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts"
SCENARIO <- "moderate"
COUNT_FILE <- file.path(BASE, "04_confidence_filter_rescuedLR_274", SCENARIO,
                        "CWOW_274_rescuedLR.confidence_filtered.counts.tsv.gz")
MAP_FILE <- file.path(BASE, "04_confidence_filter_rescuedLR_274", SCENARIO,
                      "CWOW_274_rescuedLR.transcript_to_gene.tsv.gz")
METADATA_FILE <- file.path(PROJECT, "CWOW_full_long_read_cohort_n276.csv")
OUT <- file.path(BASE, "05_DTE_DTU", "DTE_DESeq2", SCENARIO)

USE_BATCH_BOX <- FALSE  # batch_box is NOT the sequencing batch; leave FALSE
FDR <- 0.05
MIN_FL_PER_SAMPLE <- 5L
MIN_SAMPLES_WITH_FL <- 5L
MIN_GROUP_N <- 3L
RUN_COMPARISONS <- "ALL" # all 10 disease-group comparisons
SKIP_EXISTING_DTE <- TRUE # preserve completed DESeq2 results; rebuild filtering table
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

message("Reading metadata...")
meta <- fread(METADATA_FILE, check.names = FALSE)
required <- c("Study-Subject-ID", "group", "Gender", "Age", "TDP.43")
if (USE_BATCH_BOX) required <- c(required, "batch_box")
if (!all(required %in% names(meta)))
  stop("Missing metadata columns: ", paste(setdiff(required, names(meta)), collapse = ", "))
meta[, sample_id := paste0("CWOW_", sub("^CWOW_", "", trimws(as.character(`Study-Subject-ID`))))]
meta[, group := gsub("-", "_", toupper(trimws(as.character(group))))]
meta[, sex := toupper(trimws(as.character(Gender)))]
meta[sex == "MALE", sex := "XY"]
meta[sex == "FEMALE", sex := "XX"]
meta[, Age := suppressWarnings(as.numeric(Age))]
meta[, TDP43 := suppressWarnings(as.numeric(`TDP.43`))]
if (any(!is.na(meta$TDP43) & !meta$TDP43 %in% c(0, 1)))
  stop("TDP.43 must contain only 0, 1, or missing values")
if (USE_BATCH_BOX) meta[, batch_box := as.factor(batch_box)]
if (anyDuplicated(meta$sample_id)) stop("Duplicate metadata sample IDs")
if (any(!meta$sex %in% c("XX", "XY"))) stop("Unexpected sex labels")

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
if (anyNA(meta[, .(Age, sex, group)])) stop("Missing required model covariates")
if (USE_BATCH_BOX && anyNA(meta$batch_box)) stop("Missing batch_box")

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

meta_export <- meta[, .(sample_id, group, sex, Age, TDP43)]
if (USE_BATCH_BOX) meta_export[, batch_box := meta$batch_box]
fwrite(meta_export, file.path(OUT, "matched_metadata.tsv"), sep = "\t")
fwrite(as.data.table(table(meta$group, meta$sex)), file.path(OUT, "group_by_sex.tsv"), sep = "\t")

# Center age within each comparison to reduce collinearity with the intercept.
# Age remains a continuous covariate.
make_design <- function(md, outcome = "group", adjust_group = FALSE) {
  terms <- character()
  if (nlevels(md$sex) > 1L) terms <- c(terms, "sex")
  if (length(unique(md$Age_centered)) > 1L) terms <- c(terms, "Age_centered")
  if (USE_BATCH_BOX && nlevels(md$batch_box) > 1L) terms <- c(terms, "batch_box")
  if (adjust_group && nlevels(md$group) > 1L) terms <- c(terms, "group")
  terms <- c(terms, outcome)
  f <- reformulate(terms)
  mm <- model.matrix(f, data = md)
  if (qr(mm)$rank < ncol(mm))
    stop("Design is rank-deficient: ", deparse(f), ". Check group/sex/age confounding")
  if (nrow(mm) <= ncol(mm)) stop("Too few samples for design: ", deparse(f))
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
  ix <- which(meta$group %in% pair)
  md <- as.data.frame(meta[ix])
  rownames(md) <- md$sample_id
  md$group <- factor(md$group, levels = c(pair[2], pair[1]))
  md$sex <- droplevels(factor(md$sex))
  md$Age_centered <- md$Age - mean(md$Age)
  if (USE_BATCH_BOX) md$batch_box <- droplevels(factor(md$batch_box))
  group_n <- table(md$group)
  if (length(group_n) != 2L || any(group_n < MIN_GROUP_N))
    stop("Insufficient samples in ", nm, ": ", paste(names(group_n), group_n, collapse = ", "))
  f <- make_design(md)
  out <- file.path(OUT, nm)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines(deparse(f), file.path(out, "model_formula.txt"))
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
  res <- results(dds, contrast = c("group", pair[1], pair[2]), alpha = FDR)
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
  ix <- which(!is.na(meta$TDP43))
  md <- as.data.frame(meta[ix])
  rownames(md) <- md$sample_id
  md$TDP43_status <- factor(ifelse(md$TDP43 == 1, "Positive", "Negative"),
                            levels = c("Negative", "Positive"))
  md$group <- droplevels(factor(md$group))
  md$sex <- droplevels(factor(md$sex))
  md$Age_centered <- md$Age - mean(md$Age)
  if (USE_BATCH_BOX) md$batch_box <- droplevels(factor(md$batch_box))
  status_n <- table(md$TDP43_status)
  if (any(status_n < MIN_GROUP_N)) stop("Too few TDP-43 classified samples")
  f <- make_design(md, outcome = "TDP43_status", adjust_group = TRUE)
  out <- file.path(OUT, nm)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines(deparse(f), file.path(out, "model_formula.txt"))
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
