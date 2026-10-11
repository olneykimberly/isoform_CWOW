#### NEEDS TO BE EDITED!!!
# ============================================================================
# CWOW | DIRECT SQANTI3 RECLASSIFICATION OF MODERATE-CONFIDENCE ISOMT MODELS
# Run interactively in Posit Workbench (source this R script).
#
# PURPOSE
#   Reclassify the EXACT merged IsoMatch transcript models retained in the
#   moderate-confidence DTE/DTU count matrix against the reference annotation.
#   Unlike source-tracking labels, these describe the merged GTF structures.
#
# WORKFLOW
#   1. Read retained ISOMT IDs from the DTE/DTU count matrix.
#   2. Stream merged IsoMatch GTF and write transcript/exon records for those IDs.
#   3. Run SQANTI3 QC on this subset using a matching genome FASTA + GENCODE GTF.
#   4. Check classification IDs against the retained count-matrix IDs.
#   5. Write category composition table and pie/bar plots.
#
# OUTPUT DIRECTORY: OUT (set below)
#   moderate_retained_ISOMT.gtf                 Subset GTF supplied to SQANTI3
#   SQANTI3/moderate_merged_classification.txt  Direct SQANTI3 classifications
#   SQANTI3/...                                Additional SQANTI3 QC outputs
#   direct_SQANTI3_category_assignments.tsv    One row per retained ISOMT
#   direct_SQANTI3_category_summary.tsv        Counts and percentages by class
#   direct_SQANTI3_ID_audit.tsv                ID matching and coverage summary
#   Figure_01_direct_category_pie.png/.pdf     Unique ISOMT proportions
#   Figure_02_direct_category_counts.png/.pdf  Unique ISOMT counts
#
# REQUIREMENTS
#   SQANTI3 installed and runnable via a Python interpreter; genome FASTA and
#   GENCODE reference GTF must use the SAME genome assembly/contig naming as
#   the merged IsoMatch GTF. Check annotation version used upstream.
#   Set SQANTI3_SCRIPT, PYTHON, GENOME_FASTA, REFERENCE_GTF before running.
#   Set RUN_SQANTI3 <- FALSE to reuse an existing classification output.
#   No filtering or modification of DTE/DTU input counts is performed.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# --- Editable settings ------------------------------------------------------
ROOT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts'
BASE <- file.path(ROOT, 'isomatch_outputs')
MERGED_GTF <- file.path(BASE, '02_out_CWOWSamplesCohort/CWOWSamplesCohort.merged.gtf.gz')
COUNTS <- file.path(BASE, '04_confidence_filter_rescuedLR_274/moderate/CWOW_274_rescuedLR.confidence_filtered.counts.tsv.gz')
OUT <- file.path(BASE, '05_DTE_DTU/postmerge_moderate_direct_SQANTI3')
PYTHON <- '/path/to/sqanti3/environment/bin/python' # EDIT
SQANTI3_SCRIPT <- '/path/to/SQANTI3/sqanti3_qc.py'  # EDIT
GENOME_FASTA <- '/path/to/reference/GRCh38.primary_assembly.genome.fa' # EDIT
REFERENCE_GTF <- '/path/to/reference/gencode.v38.annotation.gtf' # EDIT; match upstream reference
CPUS <- 8L
CHUNK <- 150000L
DPI <- 300L
RUN_SQANTI3 <- TRUE
PREFIX <- 'moderate_merged'

for (p in c(MERGED_GTF, COUNTS)) if (!file.exists(p)) stop('Missing: ', p)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SQANTI_OUT <- file.path(OUT, 'SQANTI3')
dir.create(SQANTI_OUT, recursive=TRUE, showWarnings=FALSE)
SUBSET_GTF <- file.path(OUT, 'moderate_retained_ISOMT.gtf')
CLASS_FILE <- file.path(SQANTI_OUT, paste0(PREFIX, '_classification.txt'))
save_tsv <- function(x, filename) fwrite(x, file.path(OUT, filename), sep='\t', na='NA')
save_plot <- function(p, stem, w=11, h=7) {
  ggsave(file.path(OUT, paste0(stem, '.png')), p, width=w, height=h, dpi=DPI, bg='white')
  ggsave(file.path(OUT, paste0(stem, '.pdf')), p, width=w, height=h, bg='white')
}
COLORS <- c(FSM='#4F83B6', ISM='#E99545', NIC='#F5C342', NNC='#9D9FA1',
            Fusion='#75B1D0', Genic='#79AA5B', `Genic intron`='#B99A73',
            Antisense='#29496C', Intergenic='#963B20', `More junctions`='#C65C91',
            Other='#B8B8B8', Unclassified='#DADADA')
normalize_class <- function(x) {
  x <- tolower(trimws(as.character(x)))
  fcase(x %in% c('full-splice_match','fsm'), 'FSM',
        x %in% c('incomplete-splice_match','ism'), 'ISM',
        x %in% c('novel_in_catalog','nic'), 'NIC',
        x %in% c('novel_not_in_catalog','nnc'), 'NNC',
        grepl('fusion',x), 'Fusion',
        x %in% c('genic','genic_genomic','genic-genomic'), 'Genic',
        x %in% c('genic_intron','genic-intron'), 'Genic intron',
        grepl('antisense',x), 'Antisense',
        grepl('intergenic',x), 'Intergenic',
        grepl('more.*junction',x), 'More junctions',
        default='Other')
}

# --- 1. Read only transcript IDs from the actual moderate count matrix ------
message('Reading retained transcript IDs...')
retained <- as.character(fread(COUNTS, select=1L, showProgress=TRUE)[[1L]])
if (anyNA(retained) || anyDuplicated(retained) || !length(retained))
  stop('Missing, duplicate, or empty retained transcript IDs.')
if (any(!grepl('^ISOMT_', retained))) stop('Unexpected transcript IDs in count matrix.')
message('Retained ISOMT IDs: ', format(length(retained), big.mark=','))

# --- 2. Stream GTF and extract retained transcript and exon features --------
# Avoid reading the full merged GTF into memory.
message('Writing retained merged GTF records...')
input <- if (grepl('\\.gz$', MERGED_GTF)) gzfile(MERGED_GTF, 'rt') else file(MERGED_GTF, 'rt')
output <- file(SUBSET_GTF, 'wt')
seen_tx <- character()
n_written <- 0L
repeat {
  lines <- readLines(input, n=CHUNK, warn=FALSE)
  if (!length(lines)) break
  features <- lines[!startsWith(lines, '#')]
  features <- features[grepl('transcript_id[[:space:]]+"', features)]
  if (!length(features)) next
  ids <- sub('.*transcript_id[[:space:]]+"([^"]+)".*', '\\1', features, perl=TRUE)
  keep <- ids %chin% retained
  if (any(keep)) {
    writeLines(features[keep], output)
    n_written <- n_written + sum(keep)
    seen_tx <- c(seen_tx, unique(ids[keep]))
  }
}
close(input)
close(output)
seen_tx <- unique(seen_tx)
missing_gtf <- setdiff(retained, seen_tx)
if (length(missing_gtf)) {
  save_tsv(data.table(transcript_id=missing_gtf), 'retained_IDs_missing_from_merged_GTF.tsv')
  stop(length(missing_gtf), ' retained transcripts missing from merged GTF; aborting SQANTI3.')
}
message('Wrote ', format(n_written, big.mark=','), ' GTF features for ', length(seen_tx), ' models.')

# --- 3. Run SQANTI3 QC ------------------------------------------------------
# SQANTI3 CLI options vary by installed version. Basic positional arguments
# and -d/-o/--cpus are standard; inspect local --help if they differ.
if (RUN_SQANTI3) {
  for (p in c(PYTHON, SQANTI3_SCRIPT, GENOME_FASTA, REFERENCE_GTF))
    if (!file.exists(p)) stop('Configure a valid SQANTI3/reference path: ', p)
  args <- c(SQANTI3_SCRIPT, SUBSET_GTF, REFERENCE_GTF, GENOME_FASTA,
            '-d', SQANTI_OUT, '-o', PREFIX, '--cpus', as.character(CPUS))
  message('Running SQANTI3 QC; this may take substantial time...')
  log_file <- file.path(OUT, 'SQANTI3_run.log')
  status <- system2(PYTHON, args=args, stdout=log_file, stderr=log_file, wait=TRUE)
  if (!identical(as.integer(status), 0L))
    stop('SQANTI3 exited with status ', status, '. Check ', log_file)
}
if (!file.exists(CLASS_FILE)) {
  candidates <- list.files(SQANTI_OUT, pattern='_classification\\.txt$', recursive=TRUE, full.names=TRUE)
  stop('Expected classification not found: ', CLASS_FILE,
       '\nInspect SQANTI3 output. Candidates: ', paste(candidates, collapse='; '))
}

# --- 4. Validate classification coverage -----------------------------------
cls_header <- names(fread(CLASS_FILE, nrows=0))
if (!all(c('isoform','structural_category') %in% cls_header))
  stop('SQANTI3 classification lacks isoform or structural_category column.')
cls <- fread(CLASS_FILE, select=c('isoform','structural_category'))
setnames(cls, 'isoform', 'transcript_id')
if (anyDuplicated(cls$transcript_id)) stop('Duplicate SQANTI3 classification IDs.')
missing_cls <- setdiff(retained, cls$transcript_id)
extra_cls <- setdiff(cls$transcript_id, retained)
coverage <- data.table(metric=c('Retained ISOMT IDs', 'IDs in subset GTF',
                                  'SQANTI3 classification IDs', 'Retained IDs missing classification',
                                  'Classification IDs not retained'),
                       n=c(length(retained), length(seen_tx), nrow(cls),
                           length(missing_cls), length(extra_cls)))
save_tsv(coverage, 'direct_SQANTI3_ID_audit.tsv')
if (length(missing_cls))
  save_tsv(data.table(transcript_id=missing_cls), 'retained_IDs_missing_SQANTI3_classification.tsv')
if (length(extra_cls))
  save_tsv(data.table(transcript_id=extra_cls), 'SQANTI3_classifications_not_retained.tsv')
assignments <- merge(data.table(transcript_id=retained), cls,
                     by='transcript_id', all.x=TRUE, sort=FALSE)
assignments[, category := fifelse(is.na(structural_category), 'Unclassified',
                                 normalize_class(structural_category))]
save_tsv(assignments, 'direct_SQANTI3_category_assignments.tsv')
summary <- assignments[, .(n_transcripts=.N), by=category]
summary[, percentage := 100*n_transcripts/sum(n_transcripts)]
setorder(summary, -n_transcripts)
save_tsv(summary, 'direct_SQANTI3_category_summary.tsv')
if (length(missing_cls)) warning(length(missing_cls), ' retained transcripts lack direct SQANTI3 classification.')

# --- 5. Plot unique retained transcript models -----------------------------
plot_dt <- summary[n_transcripts>0]
plot_dt[, category := factor(category, levels=names(COLORS))]
legend_labels <- setNames(paste0(as.character(plot_dt$category), ' (',
                                 sprintf('%.1f%%', plot_dt$percentage), ')'),
                          as.character(plot_dt$category))
pie <- ggplot(plot_dt, aes(x=1, y=n_transcripts, fill=category)) +
  geom_col(width=1, color='white', linewidth=.25) + coord_polar(theta='y') +
  xlim(.5, 1.5) +
  scale_fill_manual(values=COLORS, breaks=names(legend_labels),
                    labels=legend_labels, drop=TRUE) +
  labs(title='Direct SQANTI3 classification of retained merged isoforms',
       subtitle='Unique moderate-confidence ISOMT models; one model counted once') +
  theme_void(base_size=12) + theme(plot.title=element_text(face='bold'))
save_plot(pie, 'Figure_01_direct_category_pie', 11, 7)
bar <- ggplot(plot_dt, aes(x=reorder(category, n_transcripts), y=n_transcripts, fill=category)) +
  geom_col(show.legend=FALSE) + coord_flip() +
  scale_fill_manual(values=COLORS) + scale_y_continuous(labels=label_comma()) +
  labs(title='Structural categories of retained merged transcript models',
       subtitle='Direct SQANTI3 classification; unique moderate-confidence ISOMT IDs',
       x=NULL, y='Unique transcript models') + theme_classic(base_size=12)
save_plot(bar, 'Figure_02_direct_category_counts', 11, 7)
message('DONE. Outputs: ', OUT)
message('Review direct_SQANTI3_ID_audit.tsv before using these figures.')
