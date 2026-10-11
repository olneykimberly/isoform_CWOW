# CWOW Kinnex/IsoMatch transcript retention and pre-DTE QC
# Reads compact summaries; does NOT load the 274-column count matrix.
# Requires: ggplot2, data.table, jsonlite
BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
SCENARIO <- 'moderate'
CONF_DIR <- file.path(BASE, '04_confidence_filter_rescuedLR_274')
DTE_DIR <- file.path(BASE, '05_DTE_DTU', 'DTE_DESeq2', SCENARIO)
QC_DIR <- file.path(BASE, '05_DTE_DTU', 'input_QC', SCENARIO)
OUT <- file.path(BASE, '05_DTE_DTU', 'preDTE_figures', SCENARIO)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
for (p in c('data.table','ggplot2','jsonlite')) if (!requireNamespace(p,quietly=TRUE)) stop('Install package: ',p)
library(data.table); library(ggplot2)

# Number of transcript models in the original merged GTF (from completed merge QC).
# Override if updated merge statistics become available.
N_MERGED_GTF <- 7853532L
N_UNFILTERED_COUNTS <- 7817593L

read_conf <- function(s) {
  path <- file.path(CONF_DIR, s, 'CWOW_274_rescuedLR.filter_summary.json')
  if (!file.exists(path)) return(NULL)
  x <- jsonlite::fromJSON(path)
  data.table(scenario=s, retained=as.numeric(x$retained_transcripts),
             input=as.numeric(x$input_count_rows),
             min_fl=as.numeric(x$thresholds$min_fl_reads),
             min_samples=as.numeric(x$thresholds$min_samples))
}
conf <- rbindlist(lapply(c('permissive','moderate','stringent'), read_conf), fill=TRUE)
if (!nrow(conf)) stop('No confidence filter_summary.json files found under: ',CONF_DIR)
fwrite(conf,file.path(OUT,'confidence_scenarios.tsv'),sep='\t')
if (!SCENARIO %in% conf$scenario) stop('No confidence summary for ',SCENARIO)
n_moderate <- conf[scenario==SCENARIO,retained][1]

# Finished DTE comparisons have a filter_summary.tsv written BEFORE model fitting.
# This allows plots to be updated while DESeq2 is still running.
files <- list.files(DTE_DIR, pattern='^filter_summary\\.tsv$',recursive=TRUE,full.names=TRUE)
comparison <- rbindlist(lapply(files,function(f) {
  x <- tryCatch(fread(f),error=function(e) NULL)
  if (is.null(x) || !all(c('comparison','n_transcripts_tested') %in% names(x))) return(NULL)
  x[, .(comparison, n_samples, n_transcripts_input, n_transcripts_tested)]
}),fill=TRUE)
if (nrow(comparison)) {
  comparison <- unique(comparison,by='comparison')
  fwrite(comparison,file.path(OUT,'dte_filter_comparisons.tsv'),sep='\t')
}

# Stage-level chart: merged GTF -> unfiltered count matrix -> confidence-filtered -> DTE.
stages <- data.table(stage=c('Merged GTF','Unfiltered count matrix',
                             paste0(SCENARIO,' confidence filter')),
                     n=c(N_MERGED_GTF,N_UNFILTERED_COUNTS,n_moderate))
if (nrow(comparison)) {
  # Example comparison is labelled explicitly, never implied to represent all DTE contrasts.
  ex <- if ('AD_vs_CONTROL' %in% comparison$comparison) 'AD_vs_CONTROL' else comparison$comparison[1]
  stages <- rbind(stages,data.table(stage=paste0('DTE filter: ',gsub('_',' ',ex)),
                                   n=comparison[comparison==ex,n_transcripts_tested][1]))
}
stages[, stage:=factor(stage,levels=stage)]
stages[, percent_of_unfiltered:=100*n/N_UNFILTERED_COUNTS]
fwrite(stages,file.path(OUT,'transcript_retention_stages.tsv'),sep='\t')

p1 <- ggplot(stages,aes(x=stage,y=n)) + geom_col(width=.66,fill='#326C91') +
  geom_text(aes(label=format(n,big.mark=',',scientific=FALSE)),hjust=-.08,size=3.8) +
  coord_flip(clip='off') + scale_y_continuous(labels=function(x) format(x,big.mark=',',scientific=FALSE),
                    expand=expansion(mult=c(0,.23))) +
  labs(title='Transcript retention across CWOW long-read processing',
       subtitle='Catalog filtering and statistical expression filtering are separate steps',
       x=NULL,y='Transcript models') + theme_classic(base_size=12) +
  theme(plot.margin=margin(12,45,12,12))
ggsave(file.path(OUT,'Figure_01_transcript_retention.png'),p1,width=10,height=5.2,dpi=320)
ggsave(file.path(OUT,'Figure_01_transcript_retention.pdf'),p1,width=10,height=5.2)

# Figure 02: include both original pre-filter totals as reference bars.
# The GTF model count and count-matrix row count are NOT identical universes.
# Percentages for confidence scenarios are calculated against the unfiltered
# count-matrix rows (the actual input to the confidence-filter script).
conf[, scenario := as.character(scenario)]
conf[, pct_of_count_input := 100 * retained / N_UNFILTERED_COUNTS]
scenario_levels <- c('Merged GTF', 'Unfiltered count matrix',
                     'Permissive', 'Moderate', 'Stringent')
fig2 <- rbindlist(list(
  data.table(group = c('Merged GTF', 'Unfiltered count matrix'),
             n = c(N_MERGED_GTF, N_UNFILTERED_COUNTS),
             type = 'Original',
             pct = c(NA_real_, 100)),
  conf[, .(group = tools::toTitleCase(scenario), n = retained,
           type = 'Confidence-filtered', pct = pct_of_count_input)]
), use.names = TRUE)
fig2[, group := factor(group, levels = scenario_levels)]
fig2[, label := ifelse(is.na(pct),
                       format(n, big.mark = ',', scientific = FALSE),
                       sprintf('%s\n(%.1f%% of count input)',
                               format(n, big.mark = ',', scientific = FALSE), pct))]
fwrite(fig2, file.path(OUT, 'Figure_02_plot_data.tsv'), sep = '\t')

filter_key <- paste(
  'Original merged GTF: all IsoMatch merged transcript models (before count-matrix QC).',
  'Unfiltered count matrix: original FL-support count rows, before confidence filtering.',
  'Permissive: >=2 FL reads in >=2 of 274 samples.',
  'Moderate: >=3 FL reads in >=3 of 274 samples.',
  'Stringent: >=5 FL reads in >=5 of 274 samples.',
  sep = '\n'
)
# Use actual thresholds from JSON, if present, to prevent stale captions.
key_lines <- c(
  'Original merged GTF: all IsoMatch merged transcript models (before count-matrix QC).',
  'Unfiltered count matrix: FL-support rows before cohort confidence filtering.'
)
for (s in c('permissive', 'moderate', 'stringent')) {
  z <- conf[scenario == s]
  if (nrow(z)) {
    key_lines <- c(key_lines, sprintf('%s: >=%s FL reads in >=%s samples.',
      tools::toTitleCase(s), format(z$min_fl[1], trim=TRUE),
      format(z$min_samples[1], trim=TRUE)))
  }
}
filter_key <- paste(key_lines, collapse = '\n')

p2 <- ggplot(fig2, aes(x = group, y = n, fill = type)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = label), vjust = -0.22, size = 3.1, lineheight = 0.95) +
  scale_fill_manual(values = c('Original' = '#687889',
                               'Confidence-filtered' = '#3C7B8D'),
                    name = 'Processing stage') +
  scale_y_continuous(labels = function(x) format(x, big.mark = ',', scientific = FALSE),
                     expand = expansion(mult = c(0, 0.19))) +
  labs(title = 'Figure 02. Original catalog and cohort confidence filters',
       subtitle = 'All counts are transcript models; percentages use unfiltered count-matrix rows as denominator',
       x = NULL, y = 'Number of transcript models', caption = filter_key) +
  theme_classic(base_size = 12) +
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 18, hjust = 1),
        plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.25,
                                    margin = margin(t = 12)),
        plot.margin = margin(12, 20, 16, 12))
ggsave(file.path(OUT, 'Figure_02_confidence_scenarios.png'), p2,
       width = 12, height = 7.3, dpi = 320)
ggsave(file.path(OUT, 'Figure_02_confidence_scenarios.pdf'), p2,
       width = 12, height = 7.3)

if (nrow(comparison)) {
  # Extract ordering from the table BEFORE converting the comparison column to a factor.
  comparison_levels <- comparison[order(n_transcripts_tested), as.character(comparison)]
  comparison[, comparison := factor(as.character(comparison), levels = rev(comparison_levels))]
  p3 <- ggplot(comparison,aes(x=comparison,y=n_transcripts_tested)) +
    geom_col(fill='#477D6A',width=.7) +
    geom_text(aes(label=format(n_transcripts_tested,big.mark=',',scientific=FALSE)),hjust=-.1,size=3.3) +
    coord_flip(clip='off') + scale_y_continuous(labels=function(x) format(x,big.mark=',',scientific=FALSE),
                      expand=expansion(mult=c(0,.25))) +
    labs(title='Additional expression filtering for DESeq2',
         subtitle='At least 5 FL reads in at least 5 samples for each pairwise comparison',
         x=NULL,y='Transcripts retained for testing') + theme_classic(base_size=11) +
    theme(plot.margin=margin(10,40,10,10))
  ggsave(file.path(OUT,'Figure_03_DTE_expression_filter_by_comparison.png'),p3,
         width=10,height=max(4,2.5+.38*nrow(comparison)),dpi=320)
  ggsave(file.path(OUT,'Figure_03_DTE_expression_filter_by_comparison.pdf'),p3,
         width=10,height=max(4,2.5+.38*nrow(comparison)))
}

# Sample-level QC file generated by earlier input-QC script, if present.
sample_file <- file.path(QC_DIR,'sample_count_QC.tsv')
if (file.exists(sample_file)) {
  sq <- fread(sample_file)
  id <- intersect(c('sample_id','sample'),names(sq))
  numeric_cols <- names(sq)[vapply(sq,is.numeric,logical(1))]
  total_col <- intersect(c('total_FL_counts','total_fl_counts','total_counts','total_FL','total_reads'),numeric_cols)
  if (length(id) && length(total_col)) {
    setnames(sq,c(id[1],total_col[1]),c('sample_id','total_fl'))
    p4 <- ggplot(sq,aes(x=total_fl)) + geom_histogram(bins=35,fill='#6978A8',color='white') +
      scale_x_continuous(labels=function(x) format(x,big.mark=',',scientific=FALSE)) +
      labs(title='Figure 04. Distribution of FL-read support across 274 samples',
           subtitle='Each sample contributes one observation; totals are sums of IsoMatch FL-support counts',
           x='Total aggregated FL counts per sample',y='Number of samples') +
      theme_classic(base_size=12)
    ggsave(file.path(OUT,'Figure_04_sample_FL_support.png'),p4,width=8,height=5,dpi=320)
  } else message('Sample QC found, but total FL count column name was not recognized: ',paste(names(sq),collapse=', '))
}

message('Plots and summary tables saved to: ',OUT)
message('Rerun Source as more DTE comparisons write their filter_summary.tsv files.')
