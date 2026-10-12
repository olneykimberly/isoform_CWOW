# ============================================================================
# CWOW IsoMatch/SQANTI3 FL count reconciliation -- disk-backed, resumable audit
# SOURCE IN POSIT WORKBENCH (no command-line arguments).
# READ ONLY with respect to source data and DTE files.
# Tracking is streamed in chunks and partitioned on disk by sample, avoiding
# the huge in-memory tracking grouping that terminated earlier runs.
# IMPORTANT: SQANTI3 FL may not have the same counting semantics as IsoMatch.
# A discrepancy is diagnostic, not automatically evidence of erroneous counts.
# ============================================================================

# --------------------------- USER SETTINGS ---------------------------------
BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
PREMERGE_DIR <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results'
TRACK <- file.path(BASE,'02_out_CWOWSamplesCohort/CWOWSamplesCohort.track.tsv.gz')
UNFILTERED <- file.path(BASE,'03_qc_results_CWOWSamplesCohort_rescuedLR_dedup274_v3',
                        'rescuedLR_dedup274_v3_merged_raw_count_matrix_unfiltered.tsv')
OUT <- file.path(BASE,'05_DTE_DTU','QC_IsoMatch_FL_reconciliation_disk_backed')
TRACK_LINES_PER_CHUNK <- 50000L    # reduce to 10000 if memory is constrained
N_BUCKETS <- 128L                 # larger number -> smaller mapping partitions
SAMPLES_PER_BLOCK <- 2L          # reduce to 1 for low memory
RESUME <- TRUE
TOLERANCE <- 0
MAX_EXAMPLES <- 10000L
# ---------------------------------------------------------------------------

if (!requireNamespace('data.table',quietly=TRUE)) stop('Install data.table')
if (!requireNamespace('ggplot2',quietly=TRUE)) stop('Install ggplot2')
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
LOG_FILE <- file.path(OUT,paste0('FLNC_disk_audit_',format(Sys.time(),'%Y%m%d_%H%M%S'),'.log'))
log_con <- file(LOG_FILE,'wt')
log_line <- function(...) {
  s <- paste0(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),' | ',paste0(...,collapse=''))
  cat(s,'\n'); cat(s,'\n',file=log_con);flush(log_con)
}
write_tab <- function(x,path) fwrite(x,path,sep='\t',na='NA')
map_dir <- file.path(OUT,'mapping_by_bucket')
ck_dir <- file.path(OUT,'per_sample_checkpoints')
dir.create(map_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(ck_dir,recursive=TRUE,showWarnings=FALSE)

# Sample-based deterministic bucket assignment; no large keyed tracking table.
# Mapping is partitioned by the *sample* so a source ID never spans buckets.
bucket_of <- function(sample_ids, sample_order) {
  i <- match(sample_ids,sample_order)
  if (anyNA(i)) stop('Unexpected sample ID in tracking partition')
  ((i-1L) %% N_BUCKETS)+1L
}

run_audit <- function() {
  log_line('START; log: ',LOG_FILE)
  for(p in c(TRACK,UNFILTERED)) if(!file.exists(p)) stop('Missing input: ',p)
  if(!dir.exists(PREMERGE_DIR)) stop('Missing SQANTI3 directory: ',PREMERGE_DIR)
  if(TRACK_LINES_PER_CHUNK < 1L || SAMPLES_PER_BLOCK < 1L || N_BUCKETS < 1L)
    stop('Invalid chunk/block/bucket settings')

  log_line('Indexing SQANTI3 classifications')
  paths <- list.files(PREMERGE_DIR,pattern='_rescued_classification\\.txt(\\.gz)?$',
                      recursive=TRUE,full.names=TRUE)
  paths <- paths[!grepl('with_sr',paths,ignore.case=TRUE)]
  index <- data.table(sample_id=sub('_rescued_classification\\.txt(\\.gz)?$','',basename(paths)),path=paths)
  if(!nrow(index) || anyDuplicated(index$sample_id)) stop('Missing/duplicate SQANTI3 sample IDs')
  log_line('SQANTI3 classification files: ',nrow(index))

  header <- names(fread(UNFILTERED,nrows=0L,check.names=FALSE))
  log_line('Unfiltered matrix header: ',paste(head(header,6),collapse=' | '))
  id_col <- 'merged_isoform_id'
  if(!(id_col %in% header)) stop('Expected merged_isoform_id column not found; columns: ',paste(head(header,10),collapse=', '))
  samples <- setdiff(header,id_col)
  if(length(samples)!=274L || anyDuplicated(samples)) stop('Expected 274 unique sample columns, found ',length(samples))
  if(!all(samples %in% index$sample_id)) stop('Missing SQANTI3 classifications for matrix samples: ',paste(setdiff(samples,index$sample_id),collapse=', '))
  extra <- index[!sample_id %in% samples]
  write_tab(extra,file.path(OUT,'SQANTI3_samples_not_in_count_matrix.tsv'))
  log_line('Count matrix cohort: ',length(samples),' samples; extra SQANTI3 files: ',nrow(extra),
           if(nrow(extra)) paste0(' [',paste(extra$sample_id,collapse=', '),']') else '')
  sample_order <- sort(samples)
  setkey(index,sample_id)

  # Build disk-backed mapping once. Completion flag is written only when all
  # partitions have been deduplicated successfully. Interrupted builds restart.
  marker <- file.path(map_dir,'mapping_COMPLETE.txt')
  config <- paste(normalizePath(TRACK),file.info(TRACK)$size,as.numeric(file.info(TRACK)$mtime),
                  N_BUCKETS,paste(sample_order,collapse=';'),sep='|')
  reuse <- RESUME && file.exists(marker) && identical(readLines(marker,warn=FALSE),config)
  if(!reuse) {
    log_line('Building tracking mappings on disk (no full tracking load)')
    unlink(list.files(map_dir,full.names=TRUE),recursive=TRUE)
    # Streaming compressed tracking file; retain only TRACK_LINES_PER_CHUNK rows.
    con <- if(grepl('\\.gz$',TRACK,ignore.case=TRUE)) gzfile(TRACK,'rt') else file(TRACK,'rt')
    on.exit(try(close(con),silent=TRUE),add=TRUE)
    h <- strsplit(readLines(con,n=1L,warn=FALSE),'\t',fixed=TRUE)[[1]]
    cols <- match(c('merged_tx_id','src_tx_id','src_file_name'),h)
    if(anyNA(cols)) stop('Missing required tracking columns: ',paste(h,collapse=', '))
    n_lines <- 0; n_unmatched <- 0; n_extra <- 0; chunk_id <- 0L
    repeat {
      lines <- readLines(con,n=TRACK_LINES_PER_CHUNK,warn=FALSE)
      if(!length(lines)) break
      chunk_id <- chunk_id+1L
      # fread(text=...) is bounded to one small chunk.
      z <- fread(text=paste(lines,collapse='\n'),header=FALSE,sep='\t',
                 select=cols,showProgress=FALSE,fill=FALSE,check.names=FALSE)
      if(ncol(z)!=3L) stop('Unexpected tracking chunk column count at chunk ',chunk_id)
      # fread returns columns in file order, so restore the requested names by index.
      setnames(z,h[sort(cols)])
      setcolorder(z,c('merged_tx_id','src_tx_id','src_file_name'))
      z[,sample_id:=sub('(_rescued.*|\\.sorted\\.gtf(\\.gz)?$|\\.gtf(\\.gz)?$)','',basename(as.character(src_file_name)))]
      z[,sample_id:=sub('^(CWOW_[^._/]+).*','\\1',sample_id)]
      n_lines <- n_lines+nrow(z)
      n_unmatched <- n_unmatched+sum(!z$sample_id %in% index$sample_id)
      n_extra <- n_extra+sum(z$sample_id %in% extra$sample_id)
      z <- z[sample_id %in% samples,.(sample_id,src_tx_id=as.character(src_tx_id),merged_tx_id=as.character(merged_tx_id))]
      if(nrow(z)) {
        z[,bucket:=bucket_of(sample_id,sample_order)]
        for(b in unique(z$bucket)) {
          f <- file.path(map_dir,sprintf('bucket_%03d.tsv',b))
          fwrite(z[bucket==b,.(sample_id,src_tx_id,merged_tx_id)],f,sep='\t',
                 append=file.exists(f),col.names=!file.exists(f))
        }
      }
      rm(z,lines)
      if(chunk_id %% 20L==0L) {
        log_line('Tracking streamed: ',format(n_lines,big.mark=','),' rows, ',chunk_id,' chunks; memory MB: ',
                 round(sum(gc()[,'used']*c(8,8))/1024^2,1),' (approx)')
        gc(verbose=FALSE)
      }
    }
    close(con)
    log_line('Tracking stream complete: ',n_lines,' rows; unmatched SQANTI3 rows: ',n_unmatched,
             '; rows for extra SQANTI3 samples: ',n_extra)
    if(n_lines==0L || n_unmatched/n_lines > .01) stop('Tracking-to-SQANTI3 mapping rate below 99%')

    duplicate_parts <- list();ambiguous_parts <- list();diagnostics <- list()
    for(b in seq_len(N_BUCKETS)) {
      f <- file.path(map_dir,sprintf('bucket_%03d.tsv',b))
      if(!file.exists(f)) next
      log_line('Deduplicating tracking bucket ',b,'/',N_BUCKETS)
      d <- fread(f,showProgress=FALSE,colClasses='character')
      pair_n <- d[,.(n_tracking_rows=.N,n_merged=uniqueN(merged_tx_id)),by=.(sample_id,src_tx_id)]
      duplicate_parts[[length(duplicate_parts)+1L]] <- data.table(bucket=b,
        duplicate_source_pairs=sum(pair_n$n_tracking_rows>1L),
        ambiguous_source_pairs=sum(pair_n$n_merged>1L),
        n_tracking_rows=nrow(d))
      if(any(pair_n$n_merged>1L)) {
        a <- pair_n[n_merged>1L]
        ambiguous_parts[[length(ambiguous_parts)+1L]] <- head(a,MAX_EXAMPLES)
      }
      if(any(pair_n$n_tracking_rows>1L)) {
        diagnostics[[length(diagnostics)+1L]] <- head(pair_n[n_tracking_rows>1L],
                                                      max(1L,ceiling(MAX_EXAMPLES/N_BUCKETS)))
      }
      amb <- pair_n[n_merged>1L,.(sample_id,src_tx_id)]
      d <- unique(d[,.(sample_id,src_tx_id,merged_tx_id)])
      if(nrow(amb)) d <- d[!amb,on=.(sample_id,src_tx_id)]
      # Write a separate mapping file for each sample. Since buckets are based
      # on sample ID, a sample's mappings occur in only one bucket.
      for(s in unique(d$sample_id)) {
        target <- file.path(map_dir,paste0(s,'.tsv.gz'))
        fwrite(d[sample_id==s,.(src_tx_id,merged_tx_id)],target,sep='\t')
      }
      rm(d,pair_n);gc(verbose=FALSE)
    }
    dup_stats <- rbindlist(duplicate_parts,fill=TRUE)
    write_tab(dup_stats,file.path(OUT,'tracking_bucket_summary.tsv'))
    write_tab(head(rbindlist(diagnostics,fill=TRUE),MAX_EXAMPLES),file.path(OUT,'duplicated_tracking_examples.tsv'))
    write_tab(head(rbindlist(ambiguous_parts,fill=TRUE),MAX_EXAMPLES),file.path(OUT,'ambiguous_mapping_examples.tsv'))
    total_dup <- sum(dup_stats$duplicate_source_pairs)
    total_amb <- sum(dup_stats$ambiguous_source_pairs)
    log_line('Duplicate source pairs: ',total_dup,'; ambiguous source pairs: ',total_amb)
    # Remove temporary bucket files after mapping files are complete.
    unlink(list.files(map_dir,pattern='^bucket_[0-9]+\\.tsv$',full.names=TRUE))
    writeLines(config,marker)
    log_line('Tracking mapping COMPLETE; will reuse on next run')
  } else {
    log_line('Reusing completed on-disk tracking mappings')
    dup_stats <- fread(file.path(OUT,'tracking_bucket_summary.tsv'))
    total_dup <- sum(dup_stats$duplicate_source_pairs)
    total_amb <- sum(dup_stats$ambiguous_source_pairs)
  }

  # Validate matrix IDs once, without loading 274 count columns.
  log_line('Reading merged transcript IDs only')
  ids <- as.character(fread(UNFILTERED,select=id_col,showProgress=FALSE)[[1]])
  if(anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop('Invalid or duplicate merged transcript IDs')
  log_line('Matrix transcript rows: ',length(ids))
  # A global tracking-ID overlap check is omitted to avoid reloading all mappings;
  # per-sample joins below report mapping-to-matrix overlap directly.

  blocks <- split(samples,ceiling(seq_along(samples)/SAMPLES_PER_BLOCK))
  for(b in seq_along(blocks)) {
    todo <- blocks[[b]]
    if(RESUME) todo <- todo[!vapply(todo,function(s)
      file.exists(file.path(ck_dir,paste0(s,'.tsv'))) &&
      file.exists(file.path(ck_dir,paste0(s,'_discrepancies.tsv.gz'))),logical(1))]
    if(!length(todo)) next
    log_line('Reading count matrix block ',b,'/',length(blocks),' samples: ',paste(todo,collapse=', '))
    block <- fread(UNFILTERED,select=todo,showProgress=FALSE,check.names=FALSE)
    if(nrow(block)!=length(ids)) stop('Count matrix row count changed during audit')
    for(s in todo) {
      log_line('Auditing ',s)
      p <- index[J(s),path]
      h <- names(fread(p,nrows=0L))
      fl <- intersect(c('FL','FL.BioSample_1',grep('^FL\\.',h,value=TRUE)),h)
      if(!length(fl) || !('isoform' %in% h)) stop('SQANTI3 isoform/FL column missing for ',s)
      d <- fread(p,select=c('isoform',fl[1]),showProgress=FALSE)
      setnames(d,c('isoform',fl[1]),c('src_tx_id','FL'))
      d[,`:=`(src_tx_id=as.character(src_tx_id),FL=as.numeric(FL))]
      if(anyDuplicated(d$src_tx_id) || anyNA(d$FL) || any(d$FL<0) || any(d$FL!=round(d$FL)))
        stop('Invalid SQANTI3 FL/isoform data: ',s)
      mp <- file.path(map_dir,paste0(s,'.tsv.gz'))
      t <- if(file.exists(mp)) fread(mp,colClasses='character') else
        data.table(src_tx_id=character(),merged_tx_id=character())
      if(anyDuplicated(t$src_tx_id)) stop('Source mapping not unique for ',s)
      mapped <- merge(d,t,by='src_tx_id',all.x=TRUE,sort=FALSE)
      unmapped <- mapped[is.na(merged_tx_id)]
      expected <- mapped[!is.na(merged_tx_id),.(expected_FL=sum(FL)),by=merged_tx_id]
      observed_values <- as.numeric(block[[s]])
      if(anyNA(observed_values) || any(observed_values<0) || any(observed_values!=round(observed_values)))
        stop('Invalid matrix counts for ',s)
      observed <- data.table(merged_tx_id=ids[observed_values>0],
                             observed_FL=observed_values[observed_values>0])
      # Sparse outer join: only transcripts with expected or observed nonzero FL.
      expected <- expected[expected_FL>0]
      joined <- merge(expected,observed,by='merged_tx_id',all=TRUE,sort=FALSE)
      joined[is.na(expected_FL),expected_FL:=0]
      joined[is.na(observed_FL),observed_FL:=0]
      joined[,difference:=observed_FL-expected_FL]
      mismatch <- joined[abs(difference)>TOLERANCE]
      stats <- data.table(sample_id=s,source_total_FL=sum(d$FL),
        source_unmapped_FL=sum(unmapped$FL),source_unmapped_n=nrow(unmapped),
        mapped_expected_FL=sum(expected$expected_FL),observed_total_FL=sum(observed_values),
        observed_without_expected_FL=sum(joined[expected_FL==0,observed_FL]),
        n_expected_nonzero=nrow(expected),n_observed_nonzero=nrow(observed),
        n_discrepant_transcripts=nrow(mismatch),sum_absolute_difference=sum(abs(joined$difference)),
        n_expected_ids_absent_from_matrix=sum(!expected$merged_tx_id %in% ids))
      # Completion marker (summary) written only after discrepancy file.
      fwrite(mismatch[,.(sample_id=s,merged_tx_id,expected_FL,observed_FL,difference)],
             file.path(ck_dir,paste0(s,'_discrepancies.tsv.gz')),sep='\t')
      fwrite(stats,file.path(ck_dir,paste0(s,'.tsv')),sep='\t')
      log_line('Completed ',s,'; expected=',stats$mapped_expected_FL,
               '; observed=',stats$observed_total_FL,'; discrepant IDs=',stats$n_discrepant_transcripts)
      rm(d,t,mapped,unmapped,expected,observed_values,observed,joined,mismatch,stats)
      gc(verbose=FALSE)
    }
    rm(block);gc(verbose=FALSE)
  }

  paths <- file.path(ck_dir,paste0(samples,'.tsv'))
  if(!all(file.exists(paths))) stop('Missing sample checkpoints: ',paste(samples[!file.exists(paths)],collapse=', '))
  stats <- rbindlist(lapply(paths,fread),fill=TRUE)
  write_tab(stats,file.path(OUT,'sample_reconciliation.tsv'))
  summary <- data.table(metric=c('samples','duplicate_tracking_source_pairs',
    'ambiguous_source_pairs','total_source_FL','total_unmapped_source_FL',
    'total_expected_mapped_FL','total_observed_FL','total_discrepant_sample_transcript_pairs',
    'total_absolute_count_difference','total_expected_ids_absent_from_matrix'),
    value=c(nrow(stats),total_dup,total_amb,sum(stats$source_total_FL),
            sum(stats$source_unmapped_FL),sum(stats$mapped_expected_FL),
            sum(stats$observed_total_FL),sum(stats$n_discrepant_transcripts),
            sum(stats$sum_absolute_difference),sum(stats$n_expected_ids_absent_from_matrix)))
  write_tab(summary,file.path(OUT,'overall_audit_summary.tsv'))
  p <- ggplot(stats,aes(mapped_expected_FL,observed_total_FL))+geom_point(alpha=.7)+
    geom_abline(slope=1,intercept=0,linetype=2)+theme_classic()+
    labs(title='IsoMatch FL reconciliation by sample',x='Mapped SQANTI3 FL',y='Unfiltered IsoMatch FL')
  ggsave(file.path(OUT,'Figure_01_sample_FL_reconciliation.png'),p,width=8,height=6,dpi=300)
  ggsave(file.path(OUT,'Figure_01_sample_FL_reconciliation.pdf'),p,width=8,height=6)
  p <- ggplot(stats,aes(reorder(sample_id,sum_absolute_difference),sum_absolute_difference))+
    geom_col()+coord_flip()+theme_classic()+theme(axis.text.y=element_text(size=4))+
    labs(title='Absolute FL discrepancy by sample',x='Sample',y='Absolute count difference')
  ggsave(file.path(OUT,'Figure_02_sample_discrepancy.png'),p,width=12,height=18,dpi=300)
  log_line('SUCCESS: 274-sample audit complete; inspect overall_audit_summary.tsv')
  log_line('Discrepancies do not by themselves prove incorrect IsoMatch counting semantics')
  invisible(summary)
}

tryCatch(run_audit(),error=function(e) {
  log_line('ERROR: ',conditionMessage(e))
  log_line('Audit stopped; inspect log and mapping/checkpoints')
},finally={log_line('Log saved: ',LOG_FILE);close(log_con)})
