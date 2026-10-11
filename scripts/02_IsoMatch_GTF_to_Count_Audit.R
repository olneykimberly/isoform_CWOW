# ============================================================================
# CWOW | IsoMatch merged GTF-to-unfiltered FL count matrix audit
# ============================================================================
# PURPOSE
#   Audit transcript representation between the IsoMatch merged transcript
#   catalog (GTF) and the unfiltered full-length (FL) count matrix for the
#   274-sample cohort. Identify transcript IDs missing from either file and
#   determine whether GTF-only IDs have records in IsoMatch tracking.
#   Presence/absence and tracking records do NOT establish why a transcript
#   is missing (e.g., filtering, zero FL support, or mapping issues).
#
# PIPELINE POSITION
#   Run after IsoMatch merging and unfiltered FL count matrix generation,
#   before confidence filtering and downstream DTE/DTU input QC/analysis.
#   This is a read-only audit; it does not alter upstream inputs.
#
# INPUTS (paths configured below)
#   1. CWOWSamplesCohort.merged.gtf.gz
#      IsoMatch merged GTF; only transcript features are read.
#   2. rescuedLR_dedup274_v3_merged_raw_count_matrix_unfiltered.tsv
#      Unfiltered FL count matrix; only the transcript-ID column is read.
#   3. CWOWSamplesCohort.track.tsv.gz
#      IsoMatch tracking records for merged-to-source transcript associations.
#
# OUTPUT DIRECTORY
#   05_DTE_DTU/merged_GTF_to_count_audit/ (relative to BASE below)
#
# OUTPUT TABLES (.tsv)
#   GTF_transcripts_absent_from_counts.tsv
#     Transcript IDs in the merged GTF but not in the count matrix.
#   Count_transcripts_absent_from_GTF.tsv
#     Transcript IDs in the count matrix but not in the merged GTF.
#   GTF_vs_count_summary.tsv
#     Unique transcript counts in the GTF and unfiltered count matrix.
#   Missing_transcripts_tracking_audit.tsv
#     For each GTF-only transcript: number of matching tracking rows and
#     classification as having or lacking a tracking record.
#   Missing_transcripts_tracking_categories.tsv
#     Number of GTF-only transcripts in each tracking-evidence category.
#
# OUTPUT FIGURES (each saved as .png at 300 DPI and .pdf)
#   Figure_01_GTF_vs_unfiltered_count
#     Bar chart of unique transcript counts in GTF versus count matrix.
#   Figure_02_missing_transcripts_tracking_evidence
#     Counts of missing GTF transcripts with versus without tracking records.
#   Figure_03_missing_transcripts_tracking_support  [conditional]
#     Histogram of tracking-row counts for missing GTF transcripts; produced
#     only when at least one missing transcript has a tracking record.
#
# EXECUTION / RESOURCE NOTES
#   GTF and tracking files are streamed in chunks (CHUNK below), but retaining
#   millions of transcript IDs may require several GB of memory.
#   Tracking audit currently matches the FIRST field of each tracking row;
#   confirm merged_tx_id is the first tracking column before interpretation.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})
BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
GTF <- file.path(BASE, '02_out_CWOWSamplesCohort/CWOWSamplesCohort.merged.gtf.gz')
COUNTS <- file.path(BASE, '03_qc_results_CWOWSamplesCohort_rescuedLR_dedup274_v3/rescuedLR_dedup274_v3_merged_raw_count_matrix_unfiltered.tsv')
TRACK <- file.path(BASE, '02_out_CWOWSamplesCohort/CWOWSamplesCohort.track.tsv.gz')
OUT <- file.path(BASE, '../merged_GTF_to_count_audit')
CHUNK <- 150000L
DPI <- 300
for (p in c(GTF, COUNTS, TRACK)) if (!file.exists(p)) stop('Missing file: ', p)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
write_tsv <- function(x, filename) fwrite(x, file.path(OUT, filename), sep='\t')
save_figure <- function(p, name, width=10, height=6) {
  ggsave(file.path(OUT,paste0(name,'.png')),p,width=width,height=height,dpi=DPI,bg='white')
  ggsave(file.path(OUT,paste0(name,'.pdf')),p,width=width,height=height,bg='white')
  message('Saved ',name)
}

# 1. Count-matrix transcript IDs (first column only; no 274-column matrix load)
message('Reading count-matrix transcript IDs only...')
counts_header <- names(fread(COUNTS, nrows=0, showProgress=FALSE))
message('First count-matrix column: ', counts_header[1])
count_ids <- as.character(fread(COUNTS, select=1L, showProgress=TRUE)[[1L]])
if (anyNA(count_ids) || anyDuplicated(count_ids)) stop('Count-matrix transcript IDs contain NA or duplicates')
count_ids <- unique(count_ids)
message('Count-matrix transcript IDs: ', length(count_ids))

# 2. Stream merged GTF, collecting only transcript features (not exon lines)
message('Streaming transcript features from merged GTF...')
con <- gzfile(GTF, open='rt')
parts <- list(); k <- 0L; n_gtf_features <- 0L
repeat {
  lines <- readLines(con, n=CHUNK, warn=FALSE)
  if (!length(lines)) break
  tx_lines <- lines[grepl('^[^#][^\t]*\t[^\t]*\ttranscript\t', lines)]
  if (!length(tx_lines)) next
  n_gtf_features <- n_gtf_features + length(tx_lines)
  # transcript_id is a GTF attribute, quoted in IsoMatch output
  ids <- sub('.*transcript_id[[:space:]]+"([^"]+)".*', '\\1', tx_lines, perl=TRUE)
  if (any(!grepl('^ISOMT_', ids))) stop('Unable to parse transcript_id in GTF; example: ', tx_lines[1])
  k <- k + 1L; parts[[k]] <- ids
  if (n_gtf_features %% 1000000L < length(tx_lines)) message('GTF transcript features processed: ',format(n_gtf_features,big.mark=','))
}
close(con)
gtf_ids <- unlist(parts, use.names=FALSE)
rm(parts); gc(verbose=FALSE)
if (anyDuplicated(gtf_ids)) warning('GTF has repeated transcript features; using unique transcript IDs')
gtf_ids <- unique(gtf_ids)
message('Unique GTF transcript IDs: ',length(gtf_ids))

# 3. Compare catalogs: this quantifies what is missing, not why
missing_ids <- setdiff(gtf_ids, count_ids)
count_only <- setdiff(count_ids, gtf_ids)
message('GTF IDs absent from counts: ', length(missing_ids))
message('Count IDs absent from GTF: ', length(count_only))
write_tsv(data.table(transcript_id=missing_ids),'GTF_transcripts_absent_from_counts.tsv')
write_tsv(data.table(transcript_id=count_only),'Count_transcripts_absent_from_GTF.tsv')
stages <- data.table(stage=c('Merged GTF','Unfiltered count matrix'),
                     n=c(length(gtf_ids),length(count_ids)))
write_tsv(stages,'GTF_vs_count_summary.tsv')

# 4. Stream IsoMatch tracking to see whether each missing transcript has source mappings
# This distinguishes no tracking record vs tracking record present; does not
# establish a definitive reason for absent FL counts.
message('Auditing missing transcript IDs against IsoMatch tracking...')
track_con <- gzfile(TRACK, open='rt')
track_header <- strsplit(readLines(track_con,n=1L), '\t', fixed=TRUE)[[1L]]
needed <- c('merged_tx_id','src_tx_id','src_file_name')
if (!all(needed %in% track_header)) stop('Tracking file missing columns: ',paste(setdiff(needed,track_header),collapse=', '))
ix <- match(needed,track_header)
missing_audit <- data.table(transcript_id=missing_ids, tracking_rows=0L,
                            n_distinct_source_transcripts=0L, n_distinct_source_files=0L)
# To limit memory, aggregate per chunk, then combine counts; unique source counts
# are not accumulated across chunks (source rows may repeat across chunks).
# Instead record tracking row counts and whether any source mapping exists.
missing_audit[,c('n_distinct_source_transcripts','n_distinct_source_files'):=NULL]
setkey(missing_audit, transcript_id)
processed <- 0L
repeat {
  lines <- readLines(track_con, n=CHUNK, warn=FALSE)
  if (!length(lines)) break
  # Filter lines before splitting; only ~36k target transcript IDs
  first_field <- sub('\t.*$','',lines)
  hits <- which(first_field %chin% missing_ids)
  if (length(hits)) {
    dt <- data.table(transcript_id=first_field[hits])
    chunk_summary <- dt[,.(tracking_rows=.N),by=transcript_id]
    missing_audit[chunk_summary, tracking_rows := tracking_rows + i.tracking_rows, on='transcript_id']
  }
  processed <- processed + length(lines)
  if (processed %% 3000000L < length(lines)) message('Tracking records processed: ',format(processed,big.mark=','))
}
close(track_con)
missing_audit[,tracking_evidence:=fifelse(tracking_rows>0L,
  'Has IsoMatch tracking record; no count-matrix row',
  'No IsoMatch tracking record')]
write_tsv(missing_audit,'Missing_transcripts_tracking_audit.tsv')
reason <- missing_audit[,.(n_transcripts=.N),by=tracking_evidence]
write_tsv(reason,'Missing_transcripts_tracking_categories.tsv')

# 5. Figures
p1 <- ggplot(stages, aes(x=stage,y=n,fill=stage)) +
  geom_col(width=.65,show.legend=FALSE) +
  geom_text(aes(label=comma(n)),vjust=-.4,size=4) +
  scale_y_continuous(labels=label_comma(),expand=expansion(mult=c(0,.12))) +
  scale_fill_manual(values=c('Merged GTF'='#4F83B6','Unfiltered count matrix'='#E99545')) +
  labs(title='IsoMatch merged catalog and FL-count matrix',
       subtitle=paste0(comma(length(missing_ids)),' GTF models absent from count matrix (',
                       sprintf('%.2f',100*length(missing_ids)/length(gtf_ids)),'%)'),
       x=NULL,y='Unique transcript models') + theme_classic(base_size=13)
save_figure(p1,'Figure_01_GTF_vs_unfiltered_count')

p2 <- ggplot(reason,aes(x=reorder(tracking_evidence,n_transcripts),y=n_transcripts,fill=tracking_evidence)) +
  geom_col(width=.6,show.legend=FALSE) +
  geom_text(aes(label=comma(n_transcripts)),hjust=-.1,size=4) +
  coord_flip() +
  scale_y_continuous(labels=label_comma(),expand=expansion(mult=c(0,.2))) +
  scale_fill_manual(values=c('No IsoMatch tracking record'='#963B20',
                             'Has IsoMatch tracking record; no count-matrix row'='#4F83B6')) +
  labs(title='Tracking evidence for GTF models absent from count matrix',
       subtitle='Evidence categories, not definitive removal reasons',
       x=NULL,y='Number of transcript models') + theme_classic(base_size=12)
save_figure(p2,'Figure_02_missing_transcripts_tracking_evidence',width=11,height=5)

if (any(missing_audit$tracking_rows>0L)) {
  p3 <- ggplot(missing_audit[tracking_rows>0],aes(x=tracking_rows)) +
    geom_histogram(bins=45,fill='#4F83B6',color='white') +
    scale_x_continuous(labels=label_comma()) +
    labs(title='IsoMatch tracking associations for missing transcript models',
         subtitle='Number of source tracking records per missing merged transcript',
         x='Tracking records',y='Missing merged transcript models') + theme_classic(base_size=12)
  save_figure(p3,'Figure_03_missing_transcripts_tracking_support')
}
message('DONE: ',OUT)
message('NOTE: Missing from count matrix does NOT prove SQANTI3 filtering or FL=0.')
message('To establish exact causes, audit sample-specific SQANTI3 FL and tracking joins in the count-generation pipeline.')
