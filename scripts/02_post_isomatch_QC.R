# CWOW DTE/DTU input QC 
# Outputs:
# sample_count_QC.tsv	Per-sample FL support, transcript detection, and invalid count values
# mapping_overview.tsv	Overall transcript-to-gene mapping statistics
# unmapped_transcripts.tsv	Transcripts without gene assignments (if any)
# multi_gene_transcripts.tsv	Transcripts assigned to multiple genes (if any)
# input_structure.txt	Summary of input dimensions and QC findings
# 
#

# ---- User settings ---------------------------------------------------------
BASE <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs"
SCENARIO <- "moderate" # change to "permissive" or "stringent" if desired
EXPECTED_SAMPLES <- 274L

FILTER_DIR <- file.path(BASE, "04_confidence_filter_rescuedLR_274", SCENARIO)
COUNTS_FILE <- file.path(FILTER_DIR, "CWOW_274_rescuedLR.confidence_filtered.counts.tsv.gz")
MAP_FILE <- file.path(FILTER_DIR, "CWOW_274_rescuedLR.transcript_to_gene.tsv.gz")
OUT_DIR <- file.path(BASE, "05_DTE_DTU", "input_QC", SCENARIO)

# ---- Packages --------------------------------------------------------------
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Install data.table in your Posit Workbench R environment: install.packages('data.table')")
}
library(data.table)

# ---- Validate inputs -------------------------------------------------------
for (p in c(COUNTS_FILE, MAP_FILE)) {
  if (!file.exists(p)) stop("Missing input file: ", p)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
message("Output directory: ", OUT_DIR)
message("Reading counts: ", COUNTS_FILE)
counts <- fread(COUNTS_FILE, check.names = FALSE, showProgress = TRUE)
message("Reading transcript-to-gene map: ", MAP_FILE)
map <- fread(MAP_FILE, check.names = FALSE, showProgress = TRUE)

if (ncol(counts) < 2L) stop("Count file must have transcript IDs and sample columns")
if (ncol(map) < 2L) stop("Transcript-to-gene file must have at least two columns")

# ---- Determine identifiers and sample columns -----------------------------
transcript_col <- names(counts)[1L]
tx_ids <- as.character(counts[[transcript_col]])
sample_ids <- names(counts)[-1L]
if (anyNA(tx_ids) || any(!nzchar(tx_ids))) stop("Missing transcript IDs")
if (anyDuplicated(tx_ids)) stop("Duplicate transcript IDs")
if (anyDuplicated(sample_ids)) stop("Duplicate sample names")
if (length(sample_ids) != EXPECTED_SAMPLES) {
  warning("Expected ", EXPECTED_SAMPLES, " samples; found ", length(sample_ids))
}

# ---- Per-sample QC: column-wise, avoiding a large numeric matrix copy -----
message("Computing sample QC for ", length(sample_ids), " samples...")
sample_qc <- rbindlist(lapply(seq_along(sample_ids), function(i) {
  v <- counts[[i + 1L]]
  if (!is.numeric(v)) stop("Non-numeric count column: ", sample_ids[i])
  data.table(
    sample_id = sample_ids[i],
    total_FL = sum(v, na.rm = TRUE),
    detected_transcripts = sum(v > 0, na.rm = TRUE),
    missing = sum(is.na(v)),
    negative = sum(v < 0, na.rm = TRUE),
    noninteger = sum(abs(v - round(v)) > 1e-8, na.rm = TRUE)
  )
}))
fwrite(sample_qc, file.path(OUT_DIR, "sample_count_QC.tsv"), sep = "\t")

# ---- Mapping QC ------------------------------------------------------------
# Prefer explicit column names; otherwise report schema and use first two columns.
map_names <- tolower(names(map))
tx_candidates <- c("transcript_id", "tx_id", "transcript", "tx")
gene_candidates <- c("gene_id", "gene")
tx_idx <- match(TRUE, map_names %in% tx_candidates, nomatch = 0L)
gene_idx <- match(TRUE, map_names %in% gene_candidates, nomatch = 0L)
if (tx_idx == 0L || gene_idx == 0L || tx_idx == gene_idx) {
  warning("Mapping column names not recognized; assuming first column=transcript, second=gene. Verify mapping_overview.tsv")
  tx_idx <- 1L
  gene_idx <- 2L
}
map_tx <- as.character(map[[tx_idx]])
map_gene <- as.character(map[[gene_idx]])
valid_map <- !is.na(map_tx) & nzchar(map_tx) & !is.na(map_gene) & nzchar(map_gene)
map_pairs <- unique(data.table(transcript_id = map_tx[valid_map], gene_id = map_gene[valid_map]))
map_ambiguity <- map_pairs[, .(n_genes = uniqueN(gene_id)), by = transcript_id][n_genes > 1L]
missing_tx <- tx_ids[!tx_ids %in% map_pairs$transcript_id]

mapping_overview <- data.table(
  metric = c("count_transcripts", "map_rows", "valid_unique_tx_gene_pairs",
             "count_ids_in_map", "count_ids_missing_from_map",
             "transcripts_mapped_to_multiple_genes", "map_transcript_column", "map_gene_column"),
  value = c(length(tx_ids), nrow(map), nrow(map_pairs),
            sum(tx_ids %in% map_pairs$transcript_id), length(missing_tx),
            nrow(map_ambiguity), names(map)[tx_idx], names(map)[gene_idx])
)
fwrite(mapping_overview, file.path(OUT_DIR, "mapping_overview.tsv"), sep = "\t")
if (length(missing_tx)) fwrite(data.table(transcript_id = missing_tx), file.path(OUT_DIR, "unmapped_transcripts.tsv"), sep = "\t")
if (nrow(map_ambiguity)) fwrite(map_ambiguity, file.path(OUT_DIR, "multi_gene_transcripts.tsv"), sep = "\t")

# ---- Summary ---------------------------------------------------------------
summary_lines <- c(
  paste("Scenario:", SCENARIO),
  paste("Counts file:", COUNTS_FILE),
  paste("Mapping file:", MAP_FILE),
  paste("Count dimensions:", nrow(counts), "transcripts x", length(sample_ids), "samples"),
  paste("Count transcript column:", transcript_col),
  paste("First sample columns:", paste(head(sample_ids, 12L), collapse = ", ")),
  paste("Mapping columns:", paste(names(map), collapse = ", ")),
  paste("Samples with missing counts:", sum(sample_qc$missing > 0)),
  paste("Samples with negative counts:", sum(sample_qc$negative > 0)),
  paste("Samples with noninteger counts:", sum(sample_qc$noninteger > 0)),
  paste("Transcripts absent from mapping:", length(missing_tx)),
  paste("Transcripts mapped to >1 gene:", nrow(map_ambiguity))
)
writeLines(summary_lines, file.path(OUT_DIR, "input_structure.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
message("QC complete. Files written to: ", OUT_DIR)
if (any(sample_qc$missing > 0 | sample_qc$negative > 0 | sample_qc$noninteger > 0)) {
  warning("Invalid count values detected; inspect sample_count_QC.tsv")
}
# Optionally free memory after inspection:
# rm(counts, map); gc()
