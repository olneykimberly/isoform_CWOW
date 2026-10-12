# CWOW: global FL support filter sensitivity (NO DESeq2 fitting)
# Source directly in Posit Workbench. Uses all 274 count-matrix samples.
# For each threshold pair, count transcripts with >= min FL in >= min samples.

BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
SCENARIO <- 'moderate'
COUNT_FILE <- file.path(BASE, '04_confidence_filter_rescuedLR_274', SCENARIO,
                        'CWOW_274_rescuedLR.confidence_filtered.counts.tsv.gz')
OUT <- file.path(dirname(BASE), 'DTE_SingleModelFit', SCENARIO, 'Global_Filter_Sensitivity')
FL_CUTOFFS <- c(1L, 2L, 5L, 10L, 15L, 20L, 25L)
SAMPLE_CUTOFFS <- c(1L, 2L, 5L, 10L, 15L, 20L, 25L)
REFERENCE_FL <- 5L
REFERENCE_SAMPLES <- 10L
# This script reads the counts once, then calculates all 49 combinations efficiently.
# If memory is limited, set ROW_CHUNK_SIZE lower (e.g. 10000).
ROW_CHUNK_SIZE <- 25000L

if (!requireNamespace('data.table', quietly = TRUE)) stop('Install data.table')
if (!requireNamespace('ggplot2', quietly = TRUE)) stop('Install ggplot2')
library(data.table)
library(ggplot2)
if (!file.exists(COUNT_FILE)) stop('Missing counts: ', COUNT_FILE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, paste0('filter_sensitivity_', format(Sys.time(), '%Y%m%d_%H%M%S'), '.log'))
con <- file(LOG, open = 'wt')
on.exit <- NULL
logmsg <- function(...) {
  msg <- paste0(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' | ', paste0(..., collapse = ''))
  message(msg)
  writeLines(msg, con)
  flush(con)
}
start <- Sys.time()
tryCatch({
  logmsg('Reading globally filtered count matrix: ', COUNT_FILE)
  dt <- fread(COUNT_FILE, check.names = FALSE, showProgress = TRUE)
  if (!'transcript_id' %in% names(dt)) stop('Expected transcript_id column')
  sample_cols <- setdiff(names(dt), 'transcript_id')
  if (length(sample_cols) != 274L) stop('Expected 274 samples, found ', length(sample_cols))
  if (anyDuplicated(dt$transcript_id)) stop('Duplicate transcript IDs')
  n_tx <- nrow(dt)
  n_samples <- length(sample_cols)
  logmsg('Loaded ', format(n_tx, big.mark = ','), ' transcripts x ', n_samples, ' samples')
  # For each FL threshold, tally the distribution of how many samples meet that threshold.
  # A histogram of prevalence 0..274 is sufficient for all sample-count cutoffs.
  histograms <- matrix(0, nrow = length(FL_CUTOFFS), ncol = n_samples + 1L,
                       dimnames = list(as.character(FL_CUTOFFS), as.character(0:n_samples)))
  for (from in seq.int(1L, n_tx, by = ROW_CHUNK_SIZE)) {
    to <- min(n_tx, from + ROW_CHUNK_SIZE - 1L)
    chunk <- as.matrix(dt[from:to, ..sample_cols])
    if (!is.numeric(chunk) && !is.integer(chunk)) stop('Non-numeric counts in rows ', from, '-', to)
    if (anyNA(chunk) || any(chunk < 0)) stop('Invalid counts in rows ', from, '-', to)
    for (i in seq_along(FL_CUTOFFS)) {
      prevalence <- rowSums(chunk >= FL_CUTOFFS[i])
      histograms[i, ] <- histograms[i, ] + tabulate(prevalence + 1L, nbins = n_samples + 1L)
    }
    rm(chunk)
    if (from == 1L || to == n_tx || (ceiling(to / ROW_CHUNK_SIZE) %% 10L == 0L))
      logmsg('Processed ', format(to, big.mark = ','), '/', format(n_tx, big.mark = ','), ' transcripts')
  }
  rm(dt)
  invisible(gc())
  grid <- CJ(min_fl = FL_CUTOFFS, min_samples = SAMPLE_CUTOFFS)
  grid[, n_retained := vapply(seq_len(.N), function(k) {
    sum(histograms[as.character(min_fl[k]), (min_samples[k] + 1L):(n_samples + 1L)])
  }, numeric(1))]
  grid[, n_excluded := n_tx - n_retained]
  grid[, pct_retained := 100 * n_retained / n_tx]
  grid[, pct_excluded := 100 - pct_retained]
  grid[, reference := min_fl == REFERENCE_FL & min_samples == REFERENCE_SAMPLES]
  setorder(grid, min_fl, min_samples)
  fwrite(grid, file.path(OUT, 'global_filter_sensitivity_49_combinations.tsv'), sep = '\t')
  fwrite(data.table(min_fl = rep(FL_CUTOFFS, each = n_samples + 1L),
                    n_samples_meeting_FL = rep(0:n_samples, times = length(FL_CUTOFFS)),
                    n_transcripts = as.vector(t(histograms))),
         file.path(OUT, 'transcript_prevalence_histograms.tsv'), sep = '\t')
  ref <- grid[reference == TRUE]
  logmsg('Reference filter >=', REFERENCE_FL, ' FL in >=', REFERENCE_SAMPLES, ' samples: ',
         format(ref$n_retained, big.mark = ','), '/', format(n_tx, big.mark = ','),
         ' retained (', sprintf('%.2f', ref$pct_retained), '%)')
  fwrite(data.table(total_transcripts = n_tx, total_samples = n_samples,
                    reference_min_fl = REFERENCE_FL, reference_min_samples = REFERENCE_SAMPLES,
                    reference_retained = ref$n_retained, reference_pct_retained = ref$pct_retained),
         file.path(OUT, 'filter_sensitivity_summary.tsv'), sep = '\t')

  # Heatmap: number retained; outline the current 5-FL/10-sample choice.
  plot_grid <- copy(grid)
  plot_grid[, fl_factor := factor(min_fl, levels = FL_CUTOFFS)]
  plot_grid[, sample_factor := factor(min_samples, levels = SAMPLE_CUTOFFS)]
  heat <- ggplot(plot_grid, aes(x = sample_factor, y = fl_factor, fill = pct_retained)) +
    geom_tile(color = 'white', linewidth = 0.5) +
    geom_text(aes(label = format(n_retained, big.mark = ',', scientific = FALSE),
                  color = n_retained < 200000), size = 3, show.legend = FALSE) +
    scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'white')) +
    geom_tile(data = plot_grid[reference == TRUE], fill = NA, color = 'black', linewidth = 1.6) +
    scale_fill_viridis_c(name = 'Retained (%)', option = 'viridis', limits = c(0, 100)) +
    labs(title = 'Global transcript filtering sensitivity',
         subtitle = paste0(format(n_tx, big.mark = ','), ' transcripts; ', n_samples,
                           ' samples. Labels show retained transcripts; black outline = >=5 FL in >=10 samples.'),
         x = 'Minimum number of samples', y = 'Minimum FL counts per sample') +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = 'bold'))

  curves <- ggplot(grid, aes(x = min_samples, y = n_retained, color = factor(min_fl), group = min_fl)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    geom_point(data = grid[reference == TRUE], shape = 21, size = 5, stroke = 1.2,
               color = 'black', fill = NA, inherit.aes = FALSE,
               aes(x = min_samples, y = n_retained)) +
    scale_x_continuous(breaks = SAMPLE_CUTOFFS) +
    scale_y_continuous(labels = function(x) format(x, big.mark = ',', scientific = FALSE)) +
    labs(title = 'Transcripts retained across sample-prevalence thresholds',
         subtitle = 'Each curve represents a minimum FL count per sample; circle marks the 5-FL/10-sample setting.',
         x = 'Minimum number of samples', y = 'Transcripts retained', color = 'Min FL') +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = 'bold'))

  pct <- ggplot(grid, aes(x = min_samples, y = pct_retained, color = factor(min_fl), group = min_fl)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_x_continuous(breaks = SAMPLE_CUTOFFS) +
    scale_y_continuous(labels = function(x) paste0(x, '%'), limits = c(0, 100)) +
    labs(title = 'Percentage of transcripts retained',
         x = 'Minimum number of samples', y = 'Retained transcripts (%)', color = 'Min FL') +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = 'bold'))

  plots <- list(heatmap_counts = heat, retention_curves_counts = curves,
                retention_curves_percent = pct)
  for (nm in names(plots)) {
    w <- if (nm == 'heatmap_counts') 12 else 10
    h <- if (nm == 'heatmap_counts') 7.5 else 6
    ggsave(file.path(OUT, paste0(nm, '.png')), plots[[nm]], width = w, height = h, dpi = 300, bg = 'white')
    ggsave(file.path(OUT, paste0(nm, '.pdf')), plots[[nm]], width = w, height = h, bg = 'white')
  }
  logmsg('Finished; elapsed minutes: ', round(as.numeric(difftime(Sys.time(), start, units = 'mins')), 2))
  logmsg('Results: ', OUT)
}, error = function(e) {
  logmsg('ERROR: ', conditionMessage(e))
  stop(e)
}, finally = {
  close(con)
})
