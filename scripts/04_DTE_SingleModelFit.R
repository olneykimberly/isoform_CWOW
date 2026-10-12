# Single model fit
# CWOW long-read differential transcript expression (DTE)
# Source directly in Posit Workbench; no command-line arguments.
# Memory-bounded contrast extraction in 5,000-transcript chunks.
# Resumable fitted-model checkpoint. Previous run did not save a fitted DDS,
# so its completed 11-hour fit cannot be recovered from the log alone.
# One full-cohort DESeq2 fit for all 10 disease contrasts, plus OPTIONAL
# separate TDP-43 fit restricted to classified individuals.
# IMPORTANT: Single-cohort fitting changes dispersion estimation relative to
# independently fitted pairwise models. One global prevalence filter (>=5 FL
# in >=10 of all 274 samples) defines the tested transcript universe for
# disease and TDP-43 models. BH adjustment is per contrast.

# ---------------------------- USER SETTINGS -----------------------------
BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
SCENARIO <- 'moderate'
COUNT_FILE <- file.path(BASE, '04_confidence_filter_rescuedLR_274', SCENARIO,
                        'CWOW_274_rescuedLR.confidence_filtered.counts.tsv.gz')
MAP_FILE <- file.path(BASE, '04_confidence_filter_rescuedLR_274', SCENARIO,
                      'CWOW_274_rescuedLR.transcript_to_gene.tsv.gz')
METADATA_FILE <- file.path(dirname(BASE), 'QC_technical_RIN_normalization_CCA',
                           'CWOW_n274_DTE_metadata_RIN_SMRTcell_RGID_augmented.tsv')
OUT <- file.path(dirname(BASE), 'DTE_SingleModelFit', SCENARIO)
FDR <- 0.05
MIN_FL_PER_SAMPLE <- 5L
MIN_SAMPLES_WITH_FL <- 10L
MIN_GROUP_N <- 3L
RUN_TDP43 <- TRUE  # TRUE runs an additional, separate TDP-43 model
OVERWRITE_EXISTING <- FALSE  # FALSE protects completed outputs; no unsafe reuse of fits
SAVE_FITTED_DDS <- TRUE       # ALWAYS checkpoint expensive fitted model
RESUME <- TRUE                 # Reuse validated fit checkpoint and completed comparisons
RESULTS_CHUNK_SIZE <- 5000L    # Bound peak RAM during DESeq2 contrast extraction
CHECKPOINT <- file.path(OUT, 'fitted_disease_dds.rds')
# ------------------------------------------------------------------------

for (pkg in c('data.table', 'DESeq2')) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop('Install required package: ', pkg)
}
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
for (p in c(COUNT_FILE, MAP_FILE, METADATA_FILE)) {
  if (!file.exists(p)) stop('Missing input: ', p)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUT, paste0('DTE_SingleModelFit_', format(Sys.time(), '%Y%m%d_%H%M%S'), '.log'))
log_con <- file(LOG_FILE, open = 'wt')
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' | ', paste0(..., collapse = ''))
  cat(line, '\n'); cat(line, '\n', file = log_con); flush(log_con)
}
time_step <- function(label, expr) {
  log_msg('START: ', label)
  t0 <- proc.time()[['elapsed']]
  value <- force(expr)
  log_msg('END: ', label, ' | elapsed ', round(proc.time()[['elapsed']] - t0, 1), ' sec')
  value
}
mem_report <- function(label) {
  g <- gc(verbose = FALSE)
  log_msg(label, ' | R Vcells used ', round(g['Vcells', 'used'] / 1024^2, 1),
          ' million | Vcells max ', round(g['Vcells', 'max used'] / 1024^2, 1),
          ' million | Ncells used ', round(g['Ncells', 'used'] / 1024^2, 1), ' million')
}

comparisons <- list(
  AD_vs_CONTROL = c('AD', 'CONTROL'),
  LBD_S_vs_CONTROL = c('LBD_S', 'CONTROL'),
  LBD_AS_vs_CONTROL = c('LBD_AS', 'CONTROL'),
  LBD_ATS_vs_CONTROL = c('LBD_ATS', 'CONTROL'),
  LBD_S_vs_AD = c('LBD_S', 'AD'),
  LBD_AS_vs_AD = c('LBD_AS', 'AD'),
  LBD_ATS_vs_AD = c('LBD_ATS', 'AD'),
  LBD_AS_vs_LBD_S = c('LBD_AS', 'LBD_S'),
  LBD_ATS_vs_LBD_S = c('LBD_ATS', 'LBD_S'),
  LBD_ATS_vs_LBD_AS = c('LBD_ATS', 'LBD_AS')
)

clean_pair <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(NA_character_)
  tokens <- trimws(strsplit(x, ';', fixed = TRUE)[[1]])
  tokens <- tokens[nzchar(tokens)]
  if (!length(tokens)) return(NA_character_)
  paste(sort(unique(tokens)), collapse = ';')
}
check_design <- function(md, formula) {
  mm <- model.matrix(formula, data = md)
  if (qr(mm)$rank != ncol(mm)) {
    stop('Design is rank deficient (', qr(mm)$rank, '/', ncol(mm),
         '). Inspect confounding of run_pair, group, sex, and TDP-43.')
  }
  if (nrow(mm) <= ncol(mm)) stop('No residual degrees of freedom')
  log_msg('Model matrix: ', nrow(mm), ' samples x ', ncol(mm), ' coefficients; rank ', qr(mm)$rank)
}
# Low-allocation prevalence counts: avoid creating a full transcripts x samples logical matrix.
prevalence <- function(counts, columns, threshold = MIN_FL_PER_SAMPLE) {
  n <- nrow(counts)
  present <- integer(n)
  for (j in columns) present <- present + as.integer(counts[, j] >= threshold)
  present
}
write_results <- function(dds, name, contrast, selected_ids, selected_gene_ids,
                          pair_samples, outdir, analysis) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!OVERWRITE_EXISTING && file.exists(file.path(outdir, 'DTE_results.tsv.gz'))) {
    if (RESUME && all(file.exists(file.path(outdir, c('DTE_results.tsv.gz','DTE_FDR05.tsv.gz','filter_summary.tsv','contrast_samples_metadata.tsv'))))) {
      old <- fread(file.path(outdir, 'filter_summary.tsv'))
      log_msg('SKIP completed contrast: ', name)
      return(invisible(list(n_tested=old$n_transcripts_tested[1], n_significant=old$n_significant[1])))
    }
    stop('Partial or existing output detected in ', outdir, '; inspect it before setting OVERWRITE_EXISTING=TRUE')
  }
  # For two non-reference group coefficients, DESeq2 can perform a costly
  # contrast-specific GLM operation. Checkpoint the fitted model before this.
  chunk_dir <- file.path(outdir, 'contrast_chunks')
  dir.create(chunk_dir, recursive=TRUE, showWarnings=FALSE)
  chunk_starts <- seq.int(1L, length(selected_ids), by=RESULTS_CHUNK_SIZE)
  chunk_paths <- character(length(chunk_starts))
  for (k in seq_along(chunk_starts)) {
    lo <- chunk_starts[k]
    hi <- min(lo + RESULTS_CHUNK_SIZE - 1L, length(selected_ids))
    chunk_paths[k] <- file.path(chunk_dir, sprintf('chunk_%05d_%07d_%07d.rds', k, lo, hi))
    if (RESUME && file.exists(chunk_paths[k])) {
      log_msg('RESUME existing results chunk ', k, '/', length(chunk_starts), ' for ', name)
      next
    }
    log_msg('START contrast chunk ', k, '/', length(chunk_starts), ' (',lo,'-',hi,') for ', name)
    t0 <- proc.time()[['elapsed']]
    chunk_res <- results(dds[lo:hi, , drop=FALSE], contrast=contrast,
                         alpha=FDR, independentFiltering=FALSE, cooksCutoff=FALSE)
    chunk_ans <- as.data.table(as.data.frame(chunk_res), keep.rownames='transcript_id')
    if (!identical(chunk_ans$transcript_id, selected_ids[lo:hi])) stop('Chunk IDs differ for ', name)
    tmp <- paste0(chunk_paths[k], '.writing')
    saveRDS(chunk_ans, tmp, compress=FALSE)
    if (!file.rename(tmp, chunk_paths[k])) stop('Failed to finalize chunk ', k)
    rm(chunk_res, chunk_ans)
    invisible(gc())
    log_msg('END contrast chunk ', k, '/',length(chunk_starts),' | ',
            round(proc.time()[['elapsed']]-t0,1),' sec')
  }
  ans <- time_step(paste0('Assemble contrast chunks: ',name),
                   rbindlist(lapply(chunk_paths, readRDS), use.names=TRUE))
  stopifnot(identical(ans$transcript_id, selected_ids))
  ans[, gene_id := selected_gene_ids]
  ans[, padj := p.adjust(pvalue, method = 'BH')]
  setcolorder(ans, c('transcript_id', 'gene_id', setdiff(names(ans), c('transcript_id', 'gene_id'))))
  setorder(ans, padj)
  time_step(paste0('Write output: ', name), {
    fwrite(ans, file.path(outdir, 'DTE_results.tsv.gz.writing'), sep = '\t')
    fwrite(ans[!is.na(padj) & padj < FDR], file.path(outdir, 'DTE_FDR05.tsv.gz.writing'), sep = '\t')
    fwrite(data.table(comparison = name, analysis = analysis,
                      n_samples = length(pair_samples), n_transcripts_tested = nrow(ans),
                      n_padj_nonNA = sum(!is.na(ans$padj)),
                      n_significant = sum(!is.na(ans$padj) & ans$padj < FDR),
                      min_fl = MIN_FL_PER_SAMPLE, min_samples = MIN_SAMPLES_WITH_FL),
           file.path(outdir, 'filter_summary.tsv.writing'), sep = '\t')
    for (fn in c('DTE_results.tsv.gz','DTE_FDR05.tsv.gz','filter_summary.tsv')) {
      if (!file.rename(file.path(outdir,paste0(fn,'.writing')), file.path(outdir,fn))) stop('Failed to finalize ',fn)
    }
  })
  log_msg('COMPLETED ', name, ': tested ', nrow(ans), '; FDR < ', FDR, ': ',
          sum(!is.na(ans$padj) & ans$padj < FDR))
  invisible(list(n_tested = nrow(ans), n_significant = sum(!is.na(ans$padj) & ans$padj < FDR)))
}

run_all <- function() {
  log_msg('Starting CWOW single-model DTE; log: ', LOG_FILE)
  log_msg('Resume mode: ', RESUME, '; fitted checkpoint: ', CHECKPOINT)
  if (file.exists(CHECKPOINT) && !RESUME && !OVERWRITE_EXISTING) stop('Fitted checkpoint exists; enable RESUME or explicitly overwrite')
  log_msg('Counts: ', COUNT_FILE)
  log_msg('Design: ~ 0 + group + Age_centered + RIN + run_pair + sex')
  log_msg('All 10 pairwise disease contrasts from ONE full-cohort DESeq2 fit')
  log_msg('GLOBAL filter: >=', MIN_FL_PER_SAMPLE, ' FL in >=', MIN_SAMPLES_WITH_FL, ' of ALL 274 samples; same transcript universe for disease and TDP-43')
  log_msg('TDP-43 separate model enabled: ', RUN_TDP43)
  log_msg('Parallel DESeq2 disabled to keep memory predictable; one full fit replaces 10 fits')
  log_msg('R version: ', R.version.string, '; DESeq2 ', as.character(packageVersion('DESeq2')))

  reuse_fit <- RESUME && file.exists(CHECKPOINT)
  if (reuse_fit) log_msg('Existing fitted DDS checkpoint found: skipping count import and DESeq2 fitting')
  meta <- time_step('Read and prepare metadata', {
    m <- fread(METADATA_FILE, check.names = FALSE)
    required <- c('sample_id', 'group', 'sex', 'Age', 'RIN_workbook',
                  'all_pacbio_runs', 'rgid_1', 'rgid_2', 'TDP.43')
    if (!all(required %in% names(m))) stop('Missing metadata columns: ', paste(setdiff(required, names(m)), collapse = ', '))
    m[, sample_id := trimws(as.character(sample_id))]
    m[, group := gsub('-', '_', toupper(trimws(as.character(group))))]
    m[, sex := toupper(trimws(as.character(sex)))]
    m[sex == 'MALE', sex := 'XY']; m[sex == 'FEMALE', sex := 'XX']
    m[, Age := suppressWarnings(as.numeric(Age))]
    m[, RIN := suppressWarnings(as.numeric(RIN_workbook))]
    m[, TDP43 := suppressWarnings(as.numeric(`TDP.43`))]
    if (any(!is.na(m$TDP43) & !m$TDP43 %in% c(0, 1))) stop('Invalid TDP.43 values')
    m[, run_pair_original := vapply(as.character(all_pacbio_runs), clean_pair, character(1))]
    run_keys <- sort(unique(m$run_pair_original[!is.na(m$run_pair_original)]))
    run_map <- data.table(run_pair_original = run_keys,
                          run_pair = sprintf('RunPair_%02d', seq_along(run_keys)))
    m[run_map, on = 'run_pair_original', run_pair := i.run_pair]
    fwrite(run_map, file.path(OUT, 'run_pair_safe_label_mapping.tsv'), sep = '\t')
    if (anyDuplicated(m$sample_id)) stop('Duplicate metadata IDs')
    if (any(!m$sex %in% c('XX', 'XY'))) stop('Invalid sex values')
    if (any(!m$group %in% c('CONTROL', 'AD', 'LBD_S', 'LBD_AS', 'LBD_ATS'))) stop('Invalid disease groups')
    if (anyNA(m$run_pair) || anyNA(m$rgid_1) || anyNA(m$rgid_2)) stop('Missing run pair or RGID')
    m
  })

  if (!reuse_fit) {
  raw <- time_step('Read confidence-filtered count matrix',
                   fread(COUNT_FILE, check.names = FALSE, showProgress = TRUE))
  if (names(raw)[1] != 'transcript_id') stop('First count column must be transcript_id')
  ids <- raw[['transcript_id']]
  if (anyDuplicated(ids)) stop('Duplicate transcript IDs')
  sample_names <- names(raw)[-1L]
  if (length(sample_names) != 274L) stop('Expected 274 sample columns, found ', length(sample_names))
  if (anyDuplicated(sample_names)) stop('Duplicate sample columns')
  if (!all(sample_names %in% meta$sample_id)) stop('Missing metadata for samples: ', paste(setdiff(sample_names, meta$sample_id), collapse = ', '))
  meta <- meta[match(sample_names, sample_id)]
  stopifnot(identical(meta$sample_id, sample_names))
  log_msg('Input transcripts: ', length(ids), '; samples: ', length(sample_names))

  map <- time_step('Read and align transcript-to-gene map', {
    m <- fread(MAP_FILE, select = c('transcript_id', 'gene_id'))
    if (anyDuplicated(m$transcript_id)) stop('Duplicate transcript-gene mappings')
    idx <- match(ids, m$transcript_id)
    if (anyNA(idx) || anyNA(m$gene_id[idx])) stop('Missing gene assignments')
    genes <- m$gene_id[idx]
    rm(m, idx)
    genes
  })

  counts <- time_step('Convert to integer count matrix', {
    mat <- as.matrix(raw[, -1L, with = FALSE])
    rm(raw)
    if (!is.numeric(mat) && !is.integer(mat)) stop('Count matrix is not numeric')
    if (anyNA(mat) || any(mat < 0) || any(mat != round(mat))) stop('Counts must be nonnegative integers')
    if (any(mat > .Machine$integer.max)) stop('Count exceeds integer maximum')
    storage.mode(mat) <- 'integer'
    colnames(mat) <- sample_names
    mat
  })
  # Release the original data.table before creating DESeq2 objects.
  if (exists('raw', inherits = FALSE)) rm(raw)
  invisible(gc()); mem_report('After count matrix conversion')

  keep_samples <- which(!is.na(meta$RIN) & !is.na(meta$Age))
  if (anyNA(meta$group[keep_samples]) || anyNA(meta$sex[keep_samples])) stop('Missing model covariates')
  fwrite(meta[is.na(RIN), .(sample_id, group, RIN_workbook)],
         file.path(OUT, 'samples_excluded_missing_workbook_RIN.tsv'), sep = '\t')
  log_msg('Disease model sample count: ', length(keep_samples), ' (excluded missing RIN/Age)')
  md <- as.data.frame(meta[keep_samples])
  rownames(md) <- md$sample_id
  md$group <- factor(md$group, levels = c('CONTROL', 'AD', 'LBD_S', 'LBD_AS', 'LBD_ATS'))
  md$sex <- droplevels(factor(md$sex))
  md$run_pair <- droplevels(factor(md$run_pair))
  md$Age_centered <- (md$Age - mean(md$Age)) / 10
  md$RIN <- md$RIN - mean(md$RIN)
  if (any(table(md$group) < MIN_GROUP_N)) stop('Insufficient sample size for one disease group')
  f <- ~ 0 + group + Age_centered + RIN + run_pair + sex
  check_design(md, f)
  fwrite(as.data.table(md), file.path(OUT, 'model_metadata_used.tsv'), sep = '\t')
  fwrite(data.table(group = names(table(md$group)), n = as.integer(table(md$group))),
         file.path(OUT, 'group_sample_sizes.tsv'), sep = '\t')
  writeLines(c(deparse(f), 'One full-cohort DESeq2 fit for 10 contrasts.',
               'Age_centered = (Age - cohort mean Age)/10',
               'RIN = RIN_workbook - cohort mean RIN_workbook',
               'Global filter: >=5 FL in >=10 of all 274 matrix samples, applied BEFORE fitting.',
               'Identical global transcript universe for all 10 disease contrasts and TDP-43.',
               'Benjamini-Hochberg FDR computed separately per contrast.',
               'DESeq2 independent filtering disabled; Cook cutoff disabled for consistent BH universe.'),
             file.path(OUT, 'model_formula.txt'))

  # GLOBAL prevalence filter over ALL 274 original count-matrix samples.
  # Apply BEFORE restricting to RIN-complete disease-model samples.
  n_input_transcripts <- nrow(counts)
  global_keep <- time_step('GLOBAL prevalence filter across all 274 samples',
                           prevalence(counts, seq_len(ncol(counts))) >= MIN_SAMPLES_WITH_FL)
  log_msg('Global filter retained ', sum(global_keep), ' / ', length(global_keep),
          ' transcripts (>=', MIN_FL_PER_SAMPLE, ' FL in >=', MIN_SAMPLES_WITH_FL,
          ' of 274 samples)')
  if (!any(global_keep)) stop('No transcripts pass global prevalence filter')
  selected_ids <- ids[global_keep]
  selected_gene_ids <- map[global_keep]
  fwrite(data.table(filter_scope = 'All_274_samples', n_samples_filter = ncol(counts),
                    n_transcripts_input = n_input_transcripts,
                    n_transcripts_retained = sum(global_keep),
                    min_fl_per_sample = MIN_FL_PER_SAMPLE,
                    min_samples_with_fl = MIN_SAMPLES_WITH_FL),
         file.path(OUT, 'global_transcript_filter_summary.tsv'), sep = '\t')
  # Persist the IDs once so the separate TDP-43 model uses exactly the same filter.
  fwrite(data.table(transcript_id = selected_ids, gene_id = selected_gene_ids),
         file.path(OUT, 'global_retained_transcripts.tsv.gz'), sep = '\t')
  fit_counts <- time_step('Subset globally filtered matrix to disease-model samples',
                          counts[global_keep, keep_samples, drop = FALSE])
  rownames(fit_counts) <- selected_ids
  colnames(fit_counts) <- md$sample_id
  rm(counts, map, ids, global_keep)
  invisible(gc()); mem_report('After global-filtered count matrix')

  dds <- time_step('Construct full-cohort DESeqDataSet',
                   DESeqDataSetFromMatrix(countData = fit_counts, colData = md, design = f))
  rm(fit_counts); invisible(gc()); mem_report('Before full-cohort DESeq2 fitting')
  dds <- time_step('FULL-COHORT DESeq2 FIT (one fit for all 10 contrasts)',
                   DESeq(dds, quiet = TRUE, parallel = FALSE))
  mem_report('After full-cohort DESeq2 fitting')
  log_msg('DESeq2 coefficient names: ', paste(resultsNames(dds), collapse = ', '))
  if (SAVE_FITTED_DDS) {
    # Write to a temporary file first; rename only after a complete save.
    tmp_fit <- paste0(CHECKPOINT, '.writing')
    time_step('CHECKPOINT: save full-cohort fitted DDS', saveRDS(dds, tmp_fit, compress=FALSE))
    if (!file.rename(tmp_fit, CHECKPOINT)) stop('Cannot finalize fitted DDS checkpoint')
    log_msg('CHECKPOINT READY: ', CHECKPOINT)
  }
  } else {
    dds <- time_step('RESTORE fitted disease DDS', readRDS(CHECKPOINT))
    md <- as.data.frame(colData(dds))
    if (!all(c('sample_id','group') %in% names(md))) stop('Checkpoint missing required metadata')
    if (!identical(rownames(md), colnames(dds))) stop('Checkpoint sample order mismatch')
    selected_ids <- rownames(dds)
    if (is.null(selected_ids)) stop('Checkpoint missing transcript IDs')
    selected_map <- fread(file.path(OUT, 'global_retained_transcripts.tsv.gz'))
    if (!identical(selected_ids, selected_map$transcript_id)) stop('Checkpoint transcript IDs differ from saved global filter')
    selected_gene_ids <- selected_map$gene_id
    n_input_transcripts <- fread(file.path(OUT, 'global_transcript_filter_summary.tsv'))$n_transcripts_input[1]
    sample_names <- meta$sample_id
    log_msg('RESTORED fit: ', nrow(dds), ' transcripts; ', ncol(dds), ' samples')
    log_msg('Coefficients: ', paste(resultsNames(dds), collapse=', '))
  }
  fwrite(data.table(sample_id = colnames(dds), size_factor = sizeFactors(dds)),
         file.path(OUT, 'size_factors.tsv'), sep = '\t')

  summary_rows <- list()
  for (nm in names(comparisons)) {
    pair <- comparisons[[nm]]
    coef_names <- resultsNames(dds)
    g1 <- paste0('group', pair[1]); g2 <- paste0('group', pair[2])
    if (!all(c(g1, g2) %in% coef_names)) stop('Missing coefficients: ', g1, ', ', g2)
    contrast_vec <- setNames(numeric(length(coef_names)), coef_names)
    contrast_vec[g1] <- 1; contrast_vec[g2] <- -1
    outdir <- file.path(OUT, nm)
    ix <- which(md$group %in% pair)
    dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
    fwrite(as.data.table(md[md$group %in% pair, , drop = FALSE]),
           file.path(outdir, 'contrast_samples_metadata.tsv'), sep = '\t')
    result <- write_results(dds, nm, unname(contrast_vec), selected_ids,
                            selected_gene_ids, md$sample_id[ix],
                            outdir, 'Disease_single_model')
    summary_rows[[nm]] <- data.table(analysis = 'Disease_single_model', comparison = nm,
                                      group_1 = pair[1], n_group_1 = sum(md$group == pair[1]),
                                      group_2 = pair[2], n_group_2 = sum(md$group == pair[2]),
                                      n_total = length(ix), n_transcripts_input = n_input_transcripts,
                                      n_transcripts_tested = result$n_tested,
                                      n_significant = result$n_significant,
                                      min_fl_per_sample = MIN_FL_PER_SAMPLE,
                                      min_samples_with_fl = MIN_SAMPLES_WITH_FL)
    invisible(gc())
  }
  rm(dds); invisible(gc()); mem_report('After releasing full-cohort DESeq2 object')

  if (RUN_TDP43 && RESUME && !OVERWRITE_EXISTING && all(file.exists(file.path(OUT, 'TDP43_positive_vs_negative', c('DTE_results.tsv.gz','DTE_FDR05.tsv.gz','filter_summary.tsv'))))) {
    log_msg('SKIP completed separate TDP-43 analysis')
  } else if (RUN_TDP43) {
    log_msg('START optional TDP-43 model (SECOND fit, classified samples only)')
    # TDP-43 status is not a coefficient in the full disease-only model;
    # retaining the original adjusted TDP-43 comparison requires a second fit.
    # Reread counts for the separate TDP-43 fit; use the SAME global filter
    # defined on all 274 samples, not a TDP-43-specific prevalence filter.
    tdp_ix <- which(!is.na(meta$TDP43) & !is.na(meta$RIN) & !is.na(meta$Age))
    raw_tdp <- time_step('Read only TDP-43 sample columns from counts',
                         fread(COUNT_FILE, select = c('transcript_id', meta$sample_id[tdp_ix]),
                               check.names = FALSE, showProgress = TRUE))
    tdp_md <- as.data.frame(meta[tdp_ix])
    rownames(tdp_md) <- tdp_md$sample_id
    tdp_md$group <- droplevels(factor(tdp_md$group))
    tdp_md$run_pair <- droplevels(factor(tdp_md$run_pair))
    tdp_md$sex <- droplevels(factor(tdp_md$sex))
    tdp_md$TDP43_status <- factor(ifelse(tdp_md$TDP43 == 1, 'Positive', 'Negative'),
                                  levels = c('Negative', 'Positive'))
    tdp_md$Age_centered <- (tdp_md$Age - mean(tdp_md$Age)) / 10
    tdp_md$RIN <- tdp_md$RIN - mean(tdp_md$RIN)
    if (any(table(tdp_md$TDP43_status) < MIN_GROUP_N)) stop('Too few classified TDP-43 samples')
    tf <- ~ 0 + group + Age_centered + RIN + run_pair + sex + TDP43_status
    check_design(tdp_md, tf)
    tdp_ids <- raw_tdp[['transcript_id']]
    tdp_counts <- time_step('Build TDP-43 subset matrix',
                            as.matrix(raw_tdp[, meta$sample_id[tdp_ix], with = FALSE]))
    rm(raw_tdp)
    if (anyNA(tdp_counts) || any(tdp_counts < 0) || any(tdp_counts != round(tdp_counts))) stop('Invalid TDP counts')
    storage.mode(tdp_counts) <- 'integer'
    tdp_keep <- match(selected_ids, tdp_ids)
    if (anyNA(tdp_keep)) stop('Globally retained transcript missing from TDP input')
    tdp_selected_ids <- selected_ids
    tdp_gene_ids <- selected_gene_ids
    rm(tdp_ids)
    tdp_fit_counts <- tdp_counts[tdp_keep, , drop = FALSE]
    rownames(tdp_fit_counts) <- tdp_selected_ids
    colnames(tdp_fit_counts) <- tdp_md$sample_id
    rm(tdp_counts, tdp_keep); invisible(gc())
    tdp_dds <- time_step('Construct TDP-43 DESeqDataSet',
                         DESeqDataSetFromMatrix(countData = tdp_fit_counts,
                                                colData = tdp_md, design = tf))
    rm(tdp_fit_counts); invisible(gc())
    tdp_dds <- time_step('TDP-43 DESeq2 FIT', DESeq(tdp_dds, quiet = TRUE, parallel = FALSE))
    time_step('CHECKPOINT: save TDP-43 fitted DDS', saveRDS(tdp_dds, file.path(OUT,'fitted_TDP43_dds.rds'), compress=FALSE))
    tdp_out <- file.path(OUT, 'TDP43_positive_vs_negative')
    tdp_res <- write_results(tdp_dds, 'TDP43_positive_vs_negative',
                             c('TDP43_status', 'Positive', 'Negative'),
                             tdp_selected_ids, tdp_gene_ids,
                             tdp_md$sample_id,
                             tdp_out, 'TDP43_adjusted_for_disease')
    fwrite(as.data.table(tdp_md), file.path(tdp_out, 'model_metadata_used.tsv'), sep = '\t')
    fwrite(as.data.table(table(tdp_md$group, tdp_md$TDP43_status)),
           file.path(tdp_out, 'neuropathology_by_TDP43.tsv'), sep = '\t')
    fwrite(data.table(sample_id = colnames(tdp_dds), size_factor = sizeFactors(tdp_dds)),
           file.path(tdp_out, 'size_factors.tsv'), sep = '\t')
    writeLines(deparse(tf), file.path(tdp_out, 'model_formula.txt'))
    summary_rows[['TDP43_positive_vs_negative']] <- data.table(
      analysis = 'TDP43_adjusted_for_disease', comparison = 'TDP43_positive_vs_negative',
      group_1 = 'TDP43_Positive', n_group_1 = sum(tdp_md$TDP43_status == 'Positive'),
      group_2 = 'TDP43_Negative', n_group_2 = sum(tdp_md$TDP43_status == 'Negative'),
      n_total = nrow(tdp_md), n_transcripts_input = n_input_transcripts,
      n_transcripts_tested = tdp_res$n_tested, n_significant = tdp_res$n_significant,
      min_fl_per_sample = MIN_FL_PER_SAMPLE, min_samples_with_fl = MIN_SAMPLES_WITH_FL)
    rm(tdp_dds); invisible(gc())
    log_msg('END optional TDP-43 model')
  }
  if (RUN_TDP43 && is.null(summary_rows[['TDP43_positive_vs_negative']])) {
    p <- file.path(OUT,'TDP43_positive_vs_negative','filter_summary.tsv')
    if (file.exists(p)) {
      old <- fread(p)
      summary_rows[['TDP43_positive_vs_negative']] <- data.table(
        analysis='TDP43_adjusted_for_disease', comparison='TDP43_positive_vs_negative',
        group_1='TDP43_Positive', n_group_1=NA_integer_, group_2='TDP43_Negative', n_group_2=NA_integer_,
        n_total=old$n_samples[1], n_transcripts_input=n_input_transcripts,
        n_transcripts_tested=old$n_transcripts_tested[1], n_significant=old$n_significant[1],
        min_fl_per_sample=MIN_FL_PER_SAMPLE, min_samples_with_fl=MIN_SAMPLES_WITH_FL)
    }
  }
  fwrite(rbindlist(summary_rows, use.names = TRUE, fill = TRUE),
         file.path(OUT, 'Supplementary_Table_DTE_Transcript_Filtering.tsv'), sep = '\t')
  log_msg('ALL ANALYSES COMPLETED. Results: ', OUT)
}

tryCatch(run_all(), error = function(e) {
  log_msg('ERROR: ', conditionMessage(e))
  stop(e)
}, finally = {
  log_msg('Log saved: ', LOG_FILE)
  close(log_con)
})
