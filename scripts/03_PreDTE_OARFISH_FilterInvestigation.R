# Oarfish quantification: incremental filtering sensitivity, two references.
# Source in Posit Workbench; no command-line arguments.
# No DESeq2 fits. Only finished, readable quant files are included.

# ----------------------------- USER SETTINGS -----------------------------
PREMERGE_DIR <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results'
BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
MODERATE_MAP <- file.path(BASE, '04_confidence_filter_rescuedLR_274', 'moderate',
                          'CWOW_274_rescuedLR.transcript_to_gene.tsv.gz')
OUT <- file.path(dirname(BASE), 'DTE_SingleModelFit', 'moderate', 'Oarfish_Filter_Sensitivity')
COUNT_THRESHOLDS <- c(1, 2, 5, 10, 15, 20, 25)
SAMPLE_THRESHOLDS <- c(1, 2, 5, 10, 15, 20, 25)
# Counts are estimated Oarfish read counts, not IsoMatch FL counts.
# Run on partial samples now, rerun later as more samples finish.
RESUME <- TRUE
SAVE_CHECKPOINT_EVERY <- 1L
# Set to TRUE only if you deliberately want to discard all checkpoints.
REBUILD_CHECKPOINTS <- FALSE
# -------------------------------------------------------------------------

for (p in c('data.table', 'ggplot2')) if (!requireNamespace(p, quietly = TRUE)) stop('Install ', p)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, paste0('Oarfish_filter_', format(Sys.time(), '%Y%m%d_%H%M%S'), '.log'))
con <- file(LOG, 'wt')
log_msg <- function(...) {
  x <- paste0(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' | ', paste0(..., collapse = ''))
  cat(x, '\n'); cat(x, '\n', file = con); flush(con)
}
GROUP_LEVELS <- c('Control', 'AD(AT)', 'LBD(S)', 'LBD(AS)', 'LBD(ATS)')
GROUP_COLORS <- c('Control'='#377EB8', 'AD(AT)'='#E41A1C',
                  'LBD(S)'='#6A3D9A', 'LBD(AS)'='#9E9AC8', 'LBD(ATS)'='#636363')
group_of <- function(sample) {
  suffix <- sub('^.*-', '', sample)
  # Match longest hyphenated suffix first.
  if (grepl('-CONTROL$', sample, ignore.case=TRUE)) return('Control')
  if (grepl('-LBD-ATS$', sample, ignore.case=TRUE)) return('LBD(ATS)')
  if (grepl('-LBD-AS$', sample, ignore.case=TRUE)) return('LBD(AS)')
  if (grepl('-LBD-S$', sample, ignore.case=TRUE)) return('LBD(S)')
  if (grepl('-AD-AT$', sample, ignore.case=TRUE) || grepl('-AD$', sample, ignore.case=TRUE)) return('AD(AT)')
  NA_character_
}

# Oarfish versions may differ in column labels. Require an explicit count column.
read_quant <- function(path) {
  header <- names(fread(path, nrows=0, check.names=FALSE, showProgress=FALSE))
  id_candidates <- c('name','Name','target_id','target','transcript_id','transcript','tname')
  count_candidates <- c('num_reads','NumReads','numReads','est_counts','est_count',
                        'estimated_counts','estimated_count','count','counts','read_count')
  id <- intersect(id_candidates, header)[1]
  cnt <- intersect(count_candidates, header)[1]
  if (is.na(id) || is.na(cnt)) {
    stop('Unrecognized Oarfish quantification header for ', path, ': ',
         paste(header, collapse=', '), '. Need transcript ID and estimated count columns.')
  }
  q <- fread(path, select=c(id,cnt), check.names=FALSE, showProgress=FALSE)
  setnames(q, c(id,cnt), c('transcript_id','estimated_count'))
  q[, transcript_id := as.character(transcript_id)]
  q[, estimated_count := suppressWarnings(as.numeric(estimated_count))]
  if (anyNA(q$transcript_id) || anyDuplicated(q$transcript_id)) stop('Invalid/duplicate transcript IDs: ', path)
  if (anyNA(q$estimated_count) || any(!is.finite(q$estimated_count)) || any(q$estimated_count < 0))
    stop('Invalid estimated counts: ', path)
  q
}

# Discovery accepts both *.quant files and *.quant/quant.sf directories if needed.
find_quant <- function() {
  if (!dir.exists(PREMERGE_DIR)) stop('Missing PREMERGE_DIR: ', PREMERGE_DIR)
  log_msg('Discovering completed quant files under ', PREMERGE_DIR)
  files <- list.files(PREMERGE_DIR, recursive=TRUE, full.names=TRUE,
                      pattern='(CWOWSamplesCohort_merged|gencode_v38_annotation)\\.quant$')
  files <- files[file.exists(files) & !dir.exists(files)]
  # Some quantifiers store an abundance table inside a .quant directory.
  dirs <- list.dirs(PREMERGE_DIR, recursive=TRUE, full.names=TRUE)
  dirs <- dirs[grepl('(CWOWSamplesCohort_merged|gencode_v38_annotation)\\.quant$', dirs)]
  if (length(dirs)) {
    for (d in dirs) {
      candidates <- file.path(d, c('quant.sf','quant.tsv','abundance.tsv'))
      candidates <- candidates[file.exists(candidates)]
      if (length(candidates)) files <- c(files, candidates[1])
    }
  }
  files <- unique(files)
  if (!length(files)) stop('No .quant files found. Check PREMERGE_DIR and quant file naming.')
  # Extract sample from quant basename, or enclosing directory if quant is a directory.
  nm <- basename(files)
  nested <- nm %in% c('quant.sf','quant.tsv','abundance.tsv')
  nm[nested] <- basename(dirname(files[nested]))
  kind <- ifelse(grepl('CWOWSamplesCohort_merged', nm, fixed=TRUE), 'IsoMatch_moderate', 'GENCODE_v38')
  sample <- sub('\\.(CWOWSamplesCohort_merged|gencode_v38_annotation)\\.quant$', '', nm)
  info <- file.info(files)
  tab <- data.table(path=files, sample_id=sample, reference=kind,
                    size=as.numeric(info$size), mtime=as.numeric(info$mtime))
  tab <- tab[!is.na(size) & size > 0]
  setorder(tab, reference, sample_id, -mtime)
  if (anyDuplicated(tab[, .(reference,sample_id)])) {
    log_msg('WARNING duplicate quant files for same sample/reference: using newest')
    tab <- unique(tab, by=c('reference','sample_id'))
  }
  tab[, group := vapply(sample_id, group_of, character(1))]
  fwrite(tab, file.path(OUT, 'discovered_quant_files.tsv'), sep='\t')
  tab
}

plot_outputs <- function(dt, ref_dir, reference, n_samples) {
  dt[, count_threshold := factor(count_threshold, levels=COUNT_THRESHOLDS)]
  dt[, sample_threshold := factor(sample_threshold, levels=SAMPLE_THRESHOLDS)]
  main_title <- paste0(reference, ' | ', n_samples, ' completed samples')
  p1 <- ggplot(dt, aes(x=sample_threshold, y=count_threshold, fill=n_retained)) +
    geom_tile(color='white', linewidth=0.35) +
    geom_text(aes(label=format(n_retained, big.mark=',', scientific=FALSE),
                  color=n_retained < 200000), size=3.3) +
    scale_color_manual(values=c('TRUE'='white','FALSE'='black'), guide='none') +
    scale_fill_viridis_c(option='C', labels=scales::label_comma()) +
    labs(title='Oarfish estimated-count filtering sensitivity', subtitle=main_title,
         x='Minimum number of samples', y='Minimum estimated counts per sample',
         fill='Retained transcripts') +
    theme_bw(base_size=12) + theme(panel.grid=element_blank())
  p2 <- ggplot(dt, aes(x=as.numeric(as.character(sample_threshold)), y=n_retained,
                       color=count_threshold, group=count_threshold)) +
    geom_line(linewidth=0.9) + geom_point(size=2) +
    scale_x_continuous(breaks=SAMPLE_THRESHOLDS) +
    scale_y_continuous(labels=scales::label_comma()) +
    labs(title='Transcript retention by sample prevalence', subtitle=main_title,
         x='Minimum number of samples', y='Transcripts retained',
         color='Minimum estimated counts') + theme_bw(base_size=12)
  p3 <- ggplot(dt, aes(x=as.numeric(as.character(sample_threshold)), y=pct_retained,
                       color=count_threshold, group=count_threshold)) +
    geom_line(linewidth=0.9) + geom_point(size=2) +
    scale_x_continuous(breaks=SAMPLE_THRESHOLDS) +
    labs(title='Percentage of transcripts retained', subtitle=main_title,
         x='Minimum number of samples', y='% of reference transcript universe',
         color='Minimum estimated counts') + theme_bw(base_size=12)
  plots <- list('01_Heatmap_retained_transcripts'=p1,
                '02_Retention_curves'=p2, '03_Percent_retained'=p3)
  for (nm in names(plots)) {
    ggsave(file.path(ref_dir,paste0(nm,'.png')), plots[[nm]], width=10, height=7, dpi=300)
    ggsave(file.path(ref_dir,paste0(nm,'.pdf')), plots[[nm]], width=10, height=7)
  }
}

run_reference <- function(manifest, reference, whitelist=NULL) {
  tab <- manifest[manifest$reference == reference]
  if (!nrow(tab)) { log_msg('No finished samples for ', reference); return(invisible(NULL)) }
  ref_dir <- file.path(OUT, reference)
  dir.create(ref_dir, recursive=TRUE, showWarnings=FALSE)
  cp_path <- file.path(ref_dir, 'prevalence_checkpoint.rds')
  if (REBUILD_CHECKPOINTS && file.exists(cp_path)) unlink(cp_path)
  # A whitelist is mandatory for the IsoMatch confidence-filtered analysis.
  if (identical(reference, 'IsoMatch_moderate') && is.null(whitelist)) stop('Moderate whitelist missing')
  state <- NULL
  if (RESUME && file.exists(cp_path)) {
    state <- readRDS(cp_path)
    if (!identical(state$reference,reference) || !identical(state$thresholds,COUNT_THRESHOLDS) ||
        !identical(state$whitelist_digest, if(is.null(whitelist)) NA_character_ else paste(length(whitelist),whitelist[1],tail(whitelist,1),sep='|')))
      stop('Checkpoint settings/reference changed for ', reference, '. Set REBUILD_CHECKPOINTS=TRUE.')
    log_msg('Resuming ', reference, ' with ', length(state$processed), ' processed samples')
  }
  # Establish reference ID universe from whitelist or first completed quant file.
  if (is.null(state)) {
    if (is.null(whitelist)) {
      q <- read_quant(tab$path[1]); universe <- q$transcript_id; rm(q)
    } else universe <- whitelist
    state <- list(reference=reference, thresholds=COUNT_THRESHOLDS,
                  whitelist_digest=if(is.null(whitelist)) NA_character_ else paste(length(whitelist),whitelist[1],tail(whitelist,1),sep='|'),
                  universe=universe, prevalence=matrix(0L,nrow=length(universe),ncol=length(COUNT_THRESHOLDS)),
                  processed=character(), per_sample=list(), per_sample_bins=list(), excluded=list())
    colnames(state$prevalence) <- as.character(COUNT_THRESHOLDS)
    rm(universe)
  }
  universe <- state$universe
  # Do not silently accept quant file modifications for already processed samples.
  current_keys <- paste(tab$sample_id, tab$size, tab$mtime, sep='|')
  done_samples <- sub('\\|.*$', '', state$processed)
  for (s in done_samples) {
    key <- current_keys[match(s,tab$sample_id)]
    old <- state$processed[match(s,done_samples)]
    if (is.na(key) || !identical(key,old)) stop('Previously processed quant file changed/disappeared: ', s,
                                                '. Set REBUILD_CHECKPOINTS=TRUE to recompute.')
  }
  for (i in seq_len(nrow(tab))) {
    s <- tab$sample_id[i]; key <- current_keys[i]
    if (s %in% done_samples) next
    log_msg(reference, ' | ', s, ' | sample ', i, '/', nrow(tab))
    # Failed/incomplete files are skipped, not counted as zero-expression samples.
    q <- tryCatch(read_quant(tab$path[i]), error=function(e) {
      log_msg('SKIP unreadable/incomplete ', s, ': ', conditionMessage(e)); NULL
    })
    if (is.null(q)) next
    if (is.null(whitelist)) {
      if (identical(q$transcript_id,universe)) {
        vals <- q$estimated_count
      } else {
        idx <- match(universe,q$transcript_id)
        if (anyNA(idx) || length(q$transcript_id)!=length(universe)) {
          log_msg('SKIP ', s, ': reference transcript IDs do not match established universe')
          rm(q); gc(FALSE); next
        }
        vals <- q$estimated_count[idx]
      }
    } else {
      idx <- match(universe,q$transcript_id)
      vals <- numeric(length(universe))
      found <- !is.na(idx)
      vals[found] <- q$estimated_count[idx[found]]
      # Absent transcript IDs count as zero only when the quant file is otherwise complete.
      if (mean(found) < 0.95) {
        log_msg('SKIP ', s, ': only ',round(100*mean(found),2),'% of moderate IDs present; inspect ID naming')
        rm(q,idx,vals,found); gc(FALSE); next
      }
      rm(idx,found)
    }
    for (j in seq_along(COUNT_THRESHOLDS)) {
      state$prevalence[,j] <- state$prevalence[,j] + as.integer(vals >= COUNT_THRESHOLDS[j])
    }
    state$per_sample_bins[[s]] <- data.table(sample_id=s, reference=reference,
        bin=c('0','(0,1)','[1,2)','[2,5)','[5,10)','[10,15)','[15,20)','[20,25)','25+'),
        n_transcripts=c(sum(vals==0),sum(vals>0 & vals<1),sum(vals>=1 & vals<2),
          sum(vals>=2 & vals<5),sum(vals>=5 & vals<10),sum(vals>=10 & vals<15),
          sum(vals>=15 & vals<20),sum(vals>=20 & vals<25),sum(vals>=25)))
    state$per_sample[[s]] <- data.table(sample_id=s, reference=reference,
                                        group=tab$group[i], n_reference=length(universe),
                                        n_nonzero=sum(vals>0), sum_estimated_counts=sum(vals),
                                        n_ge_1=sum(vals>=1), n_ge_5=sum(vals>=5), n_ge_10=sum(vals>=10),
                                        quant_file=tab$path[i])
    state$processed <- c(state$processed,key)
    rm(q,vals)
    if (length(state$processed) %% SAVE_CHECKPOINT_EVERY == 0L) {
      saveRDS(state, paste0(cp_path,'.tmp'), compress=FALSE)
      if (!file.rename(paste0(cp_path,'.tmp'), cp_path)) stop('Could not atomically replace checkpoint')
    }
    gc(FALSE)
  }
  saveRDS(state,paste0(cp_path,'.tmp'),compress=FALSE)
  if (!file.rename(paste0(cp_path,'.tmp'),cp_path)) stop('Could not save checkpoint')
  n_samples <- length(state$processed)
  if (!n_samples) { log_msg('No usable samples for ',reference); return(invisible(NULL)) }
  log_msg(reference, ' | usable samples: ', n_samples, '; transcript universe: ',length(universe))
  grid <- CJ(count_threshold=COUNT_THRESHOLDS, sample_threshold=SAMPLE_THRESHOLDS)
  grid[, n_retained := vapply(seq_len(.N), function(k) {
    if (sample_threshold[k] > n_samples) return(NA_integer_)
    sum(state$prevalence[,match(count_threshold[k],COUNT_THRESHOLDS)] >= sample_threshold[k])
  }, integer(1))]
  grid[, n_reference := length(universe)]
  grid[, n_samples := n_samples]
  grid[, pct_retained := 100*n_retained/n_reference]
  grid[, reference := reference]
  fwrite(grid,file.path(ref_dir,'filter_sensitivity_49_combinations.tsv'),sep='\t')
  sm <- rbindlist(state$per_sample,fill=TRUE)
  fwrite(sm,file.path(ref_dir,'per_sample_quantification_summary.tsv'),sep='\t')
  bins <- rbindlist(state$per_sample_bins,fill=TRUE)
  fwrite(bins,file.path(ref_dir,'per_sample_estimated_count_bins.tsv'),sep='\t')
  bins[, bin := factor(bin, levels=c('0','(0,1)','[1,2)','[2,5)','[5,10)','[10,15)','[15,20)','[20,25)','25+'))]
  agg <- bins[, .(mean_transcripts=mean(n_transcripts), median_transcripts=median(n_transcripts)),by=bin]
  p_dist <- ggplot(agg,aes(x=bin,y=mean_transcripts)) + geom_col(fill='#648FB4') +
    scale_y_continuous(labels=scales::label_comma()) +
    labs(title=paste0(reference,' | Estimated-count distribution'),
         subtitle=paste('Mean number of transcripts per bin across',n_samples,'samples'),
         x='Oarfish estimated counts per transcript per sample',
         y='Mean transcripts per sample') + theme_bw(base_size=12)
  ggsave(file.path(ref_dir,'05_Estimated_count_distribution.png'),p_dist,width=10,height=6,dpi=300)
  ggsave(file.path(ref_dir,'05_Estimated_count_distribution.pdf'),p_dist,width=10,height=6)
  fwrite(data.table(sample_id=tab$sample_id, included=tab$sample_id %in% sm$sample_id,
                    reference=reference,group=tab$group,path=tab$path),
         file.path(ref_dir,'sample_inclusion.tsv'),sep='\t')
  # Plot only thresholds attainable with the current number of completed samples.
  plot_dt <- copy(grid[!is.na(n_retained)])
  if (nrow(plot_dt)) plot_outputs(plot_dt,ref_dir,reference,n_samples)
  # Estimated-count sample distributions, not FL counts.
  sm[, group := factor(group,levels=GROUP_LEVELS)]
  p <- ggplot(sm[!is.na(group)],aes(x=group,y=sum_estimated_counts,fill=group)) +
    geom_boxplot(outlier.shape=NA,alpha=0.55) +
    geom_jitter(width=0.16,height=0,size=1.8) +
    scale_fill_manual(values=GROUP_COLORS,drop=FALSE) +
    scale_y_continuous(labels=scales::label_comma()) +
    labs(title=paste0(reference,' | Estimated counts per sample'),
         subtitle=paste(n_samples,'completed quantifications'),x=NULL,
         y='Sum of Oarfish estimated counts') +
    theme_bw(base_size=12) + theme(legend.position='none')
  ggsave(file.path(ref_dir,'04_Sample_total_estimated_counts_by_group.png'),p,width=9,height=6,dpi=300)
  ggsave(file.path(ref_dir,'04_Sample_total_estimated_counts_by_group.pdf'),p,width=9,height=6)
  log_msg('Finished ',reference,'; plots and tables: ',ref_dir)
  invisible(grid)
}

tryCatch({
  log_msg('Starting Oarfish filtering sensitivity; outputs: ',OUT)
  if (!file.exists(MODERATE_MAP)) stop('Moderate-confidence map not found: ',MODERATE_MAP)
  moderate <- fread(MODERATE_MAP, select='transcript_id', showProgress=TRUE)[['transcript_id']]
  if (!length(moderate) || anyDuplicated(moderate)) stop('Invalid moderate transcript whitelist')
  log_msg('Moderate-confidence transcript IDs: ',length(moderate))
  manifest <- find_quant()
  log_msg('Discovered files: ',nrow(manifest),'; merged=',sum(manifest$reference=='IsoMatch_moderate'),
          '; GENCODE=',sum(manifest$reference=='GENCODE_v38'))
  run_reference(manifest,'IsoMatch_moderate',moderate)
  rm(moderate);gc(FALSE)
  run_reference(manifest,'GENCODE_v38')
  log_msg('ALL AVAILABLE QUANTIFICATIONS ANALYZED. Rerun later to add completed samples.')
},error=function(e){log_msg('ERROR: ',conditionMessage(e));stop(e)},finally={log_msg('Log saved: ',LOG);close(con)})
