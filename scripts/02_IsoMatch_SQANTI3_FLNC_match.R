# ============================================================================
# CWOW | IsoMatch FL-support reconciliation audit | 274 samples
# PURPOSE: Reconstruct expected merged-transcript FL counts from source SQANTI3
# rescued classifications and IsoMatch tracking, then compare against the
# UNFILTERED IsoMatch count matrix. Read-only; never modifies DTE inputs.
# CAUTION: A mismatch does not automatically imply a bug: source FL support
# may differ from IsoMatch's counting inputs or source ID conventions.
# OUTPUTS: source_mapping_diagnostics.tsv, sample_reconciliation.tsv,
# transcript_discrepancies.tsv.gz, overall_audit_summary.tsv,
# unmatched_source_examples.tsv, duplicated_tracking_examples.tsv,
# ambiguous_mapping_examples.tsv, and reconciliation scatter/summary figures.
# Run directly in Posit Workbench. Adjust paths and optional limits below.
# ============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
PREMERGE_DIR <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results'
TRACK <- file.path(BASE,'02_out_CWOWSamplesCohort/CWOWSamplesCohort.track.tsv.gz')
UNFILTERED <- file.path(BASE,'03_qc_results_CWOWSamplesCohort_rescuedLR_dedup274_v3/rescuedLR_dedup274_v3_merged_raw_count_matrix_unfiltered.tsv')
OUT <- file.path(BASE,'05_DTE_DTU','QC_IsoMatch_FL_reconciliation')
MAX_EXAMPLES <- 10000L
TOLERANCE <- 0
for(p in c(TRACK,UNFILTERED)) if(!file.exists(p)) stop('Missing: ',p)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write_tab <- function(x,n) fwrite(x,file.path(OUT,n),sep='\t',na='NA')
message('Indexing original rescued SQANTI3 classification files...')
paths <- list.files(PREMERGE_DIR,pattern='_rescued_classification\\.txt(\\.gz)?$',recursive=TRUE,full.names=TRUE)
paths <- paths[!grepl('with_sr',paths,ignore.case=TRUE)]
index <- data.table(sample_id=sub('_rescued_classification\\.txt(\\.gz)?$','',basename(paths)),path=paths)
if(anyDuplicated(index$sample_id)) stop('Duplicate sample SQANTI3 files')
message('Reading tracking and checking multiplicities...')
track <- fread(TRACK,select=c('merged_tx_id','src_tx_id','src_file_name'),showProgress=TRUE)
track[,`:=`(merged_tx_id=as.character(merged_tx_id),src_tx_id=as.character(src_tx_id),src_file_name=as.character(src_file_name))]
# Source filenames may contain directory components, suffixes, or extensions.
# Verify the mapping table before interpreting discrepancies.
track[,source_basename:=basename(src_file_name)]
track[,sample_id:=sub('(_rescued.*|\\.sorted\\.gtf(\\.gz)?$|\\.gtf(\\.gz)?$)','',source_basename)]
track[,sample_id:=sub('^(CWOW_[^._/]+).*','\\1',sample_id)]
match_rate <- mean(track$sample_id %in% index$sample_id)
message('Tracking source filename -> SQANTI3 sample match fraction: ',signif(match_rate,4))
if(match_rate < .99) {
  write_tab(unique(track[!sample_id %in% index$sample_id,.(src_file_name,sample_id)])[1:min(.N,MAX_EXAMPLES)],'unmatched_tracking_source_filenames.tsv')
  stop('Tracking source filename parsing is not sufficiently reliable (',round(100*match_rate,2),'%). Inspect unmatched_tracking_source_filenames.tsv and edit sample_id parsing; do not infer count validity.')
}
setkey(track,sample_id,src_tx_id)
dup <- track[,.(n_tracking_rows=.N,n_merged=uniqueN(merged_tx_id)),by=.(sample_id,src_tx_id)]
write_tab(head(dup[n_tracking_rows>1L],MAX_EXAMPLES),'duplicated_tracking_examples.tsv')
write_tab(head(dup[n_merged>1L],MAX_EXAMPLES),'ambiguous_mapping_examples.tsv')
# Duplicate identical rows can be deduplicated safely. A source mapping to
# MULTIPLE merged IDs is not safely resolvable without IsoMatch semantics.
track <- unique(track[,.(sample_id,src_tx_id,merged_tx_id)])
amb <- track[,.(n_merged=uniqueN(merged_tx_id)),by=.(sample_id,src_tx_id)][n_merged>1L]
if(nrow(amb)) message('WARNING: ',nrow(amb),' source transcripts map to >1 merged ID. Excluding from expected count reconstruction; report as unresolved.')
track <- track[!amb,on=.(sample_id,src_tx_id)]
message('Inspecting unfiltered count matrix header...')
header <- names(fread(UNFILTERED, nrows = 0L, check.names = FALSE))
message('  First columns: ', paste(head(header, 12L), collapse = ' | '))
# Select the ID column by its name, or by validating its values against
# IsoMatch merged transcript IDs in tracking. Never assume the first column
# is transcript_id, and never accidentally use a sample-count column as IDs.
sample_cols <- intersect(header, index$sample_id)
if (length(sample_cols) != nrow(index) ||
    !setequal(sample_cols, index$sample_id)) {
  stop('Count matrix sample columns do not match the SQANTI3 sample index. ',
       'Found ', length(sample_cols), ' of ', nrow(index), ' samples. ',
       'Missing: ', paste(head(setdiff(index$sample_id, sample_cols), 8L), collapse=', '),
       '; unexpected columns: ',
       paste(head(setdiff(header, c(sample_cols, 'transcript_id')), 8L), collapse=', '))
}
id_candidates <- setdiff(header, sample_cols)
if (!length(id_candidates)) stop('No potential transcript ID column in count matrix')
known_names <- c('transcript_id','merged_tx_id','isoform','transcript','transcript_name',
                 'target_id','feature_id','ID','id','V1')
ordered_candidates <- c(intersect(known_names, id_candidates),
                        setdiff(id_candidates, known_names))
# A small diagnostic is exported regardless of the outcome.
track_ids <- unique(track$merged_tx_id)
diag <- rbindlist(lapply(ordered_candidates, function(nm) {
  probe <- fread(UNFILTERED, select = nm, showProgress = FALSE)
  vals <- as.character(probe[[1L]])
  data.table(column = nm, n_rows = length(vals),
             n_unique = uniqueN(vals), n_missing = sum(is.na(vals) | !nzchar(vals)),
             n_matching_tracking_ids = sum(vals %in% track_ids),
             fraction_matching_tracking_ids = mean(vals %in% track_ids))
}))
write_tab(diag, 'count_matrix_identifier_diagnostics.tsv')
setorder(diag, -fraction_matching_tracking_ids, -n_matching_tracking_ids)
print(diag)
if (!nrow(diag) || diag$fraction_matching_tracking_ids[1L] < 0.5) {
  stop('Could not validate an unfiltered count-matrix transcript ID column ',
       'against tracking merged IDs. Review count_matrix_identifier_diagnostics.tsv. ',
       'No reconciliation was performed.')
}
if (nrow(diag) > 1L &&
    diag$fraction_matching_tracking_ids[2L] == diag$fraction_matching_tracking_ids[1L]) {
  stop('Ambiguous count-matrix ID columns; review count_matrix_identifier_diagnostics.tsv')
}
id_col <- diag$column[1L]
message('Validated count matrix ID column: ', id_col,
        ' (', round(100*diag$fraction_matching_tracking_ids[1L], 2),
        '% match tracking IDs)')
message('Reading unfiltered count matrix (may require substantial memory)...')
raw <- fread(UNFILTERED, select = c(id_col, sample_cols), showProgress = TRUE,
             check.names = FALSE)
setnames(raw, id_col, 'transcript_id')
raw[, transcript_id := as.character(transcript_id)]
if (anyNA(raw$transcript_id) || any(!nzchar(raw$transcript_id)))
  stop('Missing or blank merged transcript IDs in count matrix')
if (anyDuplicated(raw$transcript_id))
  stop('Duplicate merged transcript IDs in unfiltered count matrix')
samples <- names(raw)[-1L]
if (!setequal(samples, index$sample_id)) stop('SQANTI3 files do not match count columns')
setkey(index,sample_id)
# Construct long expected values one sample at a time; compare only nonzero
# expected values, and additionally measure observed counts with no mapping.
result_parts <- vector('list',length(samples)); stat_parts <- vector('list',length(samples))
for(i in seq_along(samples)) {
  s <- samples[i]; p <- index[J(s),path]
  h <- names(fread(p,nrows=0))
  fl <- intersect(c('FL','FL.BioSample_1',grep('^FL\\.',h,value=TRUE)),h)
  if(!length(fl)) stop('No SQANTI3 FL column for ',s,' in ',p)
  d <- fread(p,select=c('isoform',fl[1]),showProgress=FALSE)
  setnames(d,fl[1],'FL')
  if(anyDuplicated(d$isoform)) stop('Duplicate SQANTI3 isoform IDs for ',s)
  d[,`:=`(src_tx_id=as.character(isoform),FL=as.numeric(FL))]
  if(anyNA(d$FL)|any(d$FL<0)|any(d$FL!=round(d$FL))) stop('Invalid source FL support for ',s)
  src_total <- sum(d$FL)
  t <- track[J(s),.(src_tx_id,merged_tx_id),nomatch=0L]
  mapped <- merge(d[,.(src_tx_id,FL)],t,by='src_tx_id',all.x=TRUE,sort=FALSE)
  unmapped <- mapped[is.na(merged_tx_id)]
  expected <- mapped[!is.na(merged_tx_id),.(expected_FL=sum(FL)),by=merged_tx_id]
  observed <- data.table(merged_tx_id=as.character(raw$transcript_id),observed_FL=as.numeric(raw[[s]]))
  if(anyNA(observed$observed_FL)|any(observed$observed_FL<0)) stop('Invalid count matrix for ',s)
  joined <- merge(expected,observed,by='merged_tx_id',all=TRUE,sort=FALSE)
  joined[is.na(expected_FL),expected_FL:=0]
  joined[is.na(observed_FL),observed_FL:=0]
  joined[,difference:=observed_FL-expected_FL]
  mismatch <- joined[abs(difference)>TOLERANCE & (expected_FL>0 | observed_FL>0)]
  stat_parts[[i]] <- data.table(sample_id=s,source_total_FL=src_total,
    source_unmapped_FL=sum(unmapped$FL),source_unmapped_n=nrow(unmapped),
    mapped_expected_FL=sum(expected$expected_FL),observed_total_FL=sum(observed$observed_FL),
    observed_without_expected_FL=sum(joined[expected_FL==0]$observed_FL),
    n_expected_nonzero=sum(joined$expected_FL>0),n_observed_nonzero=sum(joined$observed_FL>0),
    n_discrepant_transcripts=nrow(mismatch),sum_absolute_difference=sum(abs(joined$difference)))
  result_parts[[i]] <- mismatch[,.(sample_id=s,merged_tx_id,expected_FL,observed_FL,difference)]
  if(i%%10L==0L || i==length(samples)) message('Audited ',i,' / ',length(samples),' samples')
}
stats <- rbindlist(stat_parts); discrepancies <- rbindlist(result_parts)
write_tab(stats,'sample_reconciliation.tsv')
fwrite(discrepancies,file.path(OUT,'transcript_discrepancies.tsv.gz'),sep='\t')
summary <- data.table(metric=c('samples','source_transcripts_ambiguous_mapping','duplicate_tracking_source_pairs','total_source_FL','total_unmapped_source_FL','total_expected_mapped_FL','total_observed_FL','total_discrepant_sample_transcript_pairs','total_absolute_count_difference'),
 value=c(length(samples),nrow(amb),sum(dup$n_tracking_rows>1),sum(stats$source_total_FL),sum(stats$source_unmapped_FL),sum(stats$mapped_expected_FL),sum(stats$observed_total_FL),nrow(discrepancies),sum(stats$sum_absolute_difference)))
write_tab(summary,'overall_audit_summary.tsv')
p <- ggplot(stats,aes(mapped_expected_FL,observed_total_FL))+geom_point(alpha=.7)+geom_abline(slope=1,intercept=0,linetype=2)+theme_classic()+labs(title='IsoMatch FL-support reconciliation by sample',x='Reconstructed FL support (mapped sources)',y='Unfiltered IsoMatch FL counts')
ggsave(file.path(OUT,'Figure_01_sample_FL_reconciliation.png'),p,width=8,height=6,dpi=300)
ggsave(file.path(OUT,'Figure_01_sample_FL_reconciliation.pdf'),p,width=8,height=6)
p <- ggplot(stats,aes(reorder(sample_id,sum_absolute_difference),sum_absolute_difference))+geom_col()+coord_flip()+theme_classic()+theme(axis.text.y=element_text(size=4))+labs(title='Absolute count discrepancy by sample',x='Sample',y='Sum of absolute differences')
ggsave(file.path(OUT,'Figure_02_sample_discrepancy.png'),p,width=12,height=18,dpi=300)
message('Audit complete. Examine tracking multiplicities, source FL semantics, and mismatch tables before interpreting any discrepancy.')
