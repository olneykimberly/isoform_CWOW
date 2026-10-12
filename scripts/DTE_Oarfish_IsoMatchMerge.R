# CWOW Oarfish IsoMatch-merged moderate-confidence DTE
# Source directly in Posit Workbench. No command-line arguments.
# Quantification is against the full IsoMatch reference; filter to moderate
# transcript IDs first, then >=5 estimated counts in >=10 available samples.
# DESeq2 requires integer counts: round Oarfish estimated counts once, AFTER
# prevalence filtering on unrounded estimates. Keep this approximation explicit.
# One full-cohort fit for 10 disease contrasts; TDP-43 is a SEPARATE fit.
# IMPORTANT: partial-cohort results are preliminary; rerun when cohort completes.

# ---------------------------- EDITABLE SETTINGS --------------------------
PREMERGE_DIR <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results'
BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
MAP_FILE <- file.path(BASE, '04_confidence_filter_rescuedLR_274', 'moderate',
                      'CWOW_274_rescuedLR.transcript_to_gene.tsv.gz')
METADATA_FILE <- file.path(dirname(BASE), 'QC_technical_RIN_normalization_CCA',
                           'CWOW_n274_DTE_metadata_RIN_SMRTcell_RGID_augmented.tsv')
OUT_ROOT <- file.path(dirname(BASE), 'DTE_Oarfish_IsoMatchMerge')
MIN_EST_COUNTS <- 5
MIN_SAMPLES <- 10L
FDR <- 0.05
MIN_GROUP_N <- 3L
RUN_TDP43 <- TRUE
RESULTS_CHUNK_SIZE <- 5000L
# Set TRUE only when you intentionally want a new snapshot from updated .quant files.
# Each run writes to its OWN immutable snapshot folder (does not mix cohorts).
SNAPSHOT_LABEL <- format(Sys.time(), '%Y%m%d_%H%M%S')
# -------------------------------------------------------------------------

for (p in c('data.table', 'DESeq2'))
  if (!requireNamespace(p, quietly=TRUE)) stop('Install required package: ',p)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
for (p in c(PREMERGE_DIR, MAP_FILE, METADATA_FILE))
  if (!file.exists(p)) stop('Missing input: ',p)
OUT <- file.path(OUT_ROOT, paste0('snapshot_', SNAPSHOT_LABEL))
if (dir.exists(OUT)) stop('Snapshot output already exists: ',OUT)
dir.create(OUT, recursive=TRUE)
LOG_FILE <- file.path(OUT, 'DTE_Oarfish.log')
LOG_CON <- file(LOG_FILE,'wt')
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),' | ',paste0(...,collapse=''))
  cat(msg,'\n'); cat(msg,'\n',file=LOG_CON); flush(LOG_CON)
}
timed <- function(label, expr) {
  log_msg('START ',label); start <- proc.time()[['elapsed']]
  ans <- force(expr)
  log_msg('END ',label,' | ',round(proc.time()[['elapsed']]-start,1),' seconds')
  ans
}
mem <- function(label) {
  x <- gc()
  log_msg(label, ' | Vcells used=',round(x['Vcells','used']/1024^2,1),' million')
}
GROUPS <- c('CONTROL','AD','LBD_S','LBD_AS','LBD_ATS')
GROUP_LABELS <- c(CONTROL='Control',AD='AD(AT)',LBD_S='LBD(S)',
                  LBD_AS='LBD(AS)',LBD_ATS='LBD(ATS)')
comparisons <- list(
  AD_vs_CONTROL=c('AD','CONTROL'),
  LBD_S_vs_CONTROL=c('LBD_S','CONTROL'),
  LBD_AS_vs_CONTROL=c('LBD_AS','CONTROL'),
  LBD_ATS_vs_CONTROL=c('LBD_ATS','CONTROL'),
  LBD_S_vs_AD=c('LBD_S','AD'),
  LBD_AS_vs_AD=c('LBD_AS','AD'),
  LBD_ATS_vs_AD=c('LBD_ATS','AD'),
  LBD_AS_vs_LBD_S=c('LBD_AS','LBD_S'),
  LBD_ATS_vs_LBD_S=c('LBD_ATS','LBD_S'),
  LBD_ATS_vs_LBD_AS=c('LBD_ATS','LBD_AS'))
clean_pair <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(NA_character_)
  v <- trimws(strsplit(x,';',fixed=TRUE)[[1]])
  v <- v[nzchar(v)]
  if (!length(v)) return(NA_character_)
  paste(sort(unique(v)),collapse=';')
}
read_quant <- function(path) {
  h <- names(fread(path,nrows=0,check.names=FALSE,showProgress=FALSE))
  id <- intersect(c('name','Name','target_id','target','transcript_id','transcript','tname'),h)[1]
  ct <- intersect(c('num_reads','NumReads','numReads','est_counts','est_count',
                    'estimated_counts','estimated_count','count','counts','read_count'),h)[1]
  if (is.na(id) || is.na(ct)) stop('Unexpected quant header in ',path,': ',paste(h,collapse=', '))
  q <- fread(path,select=c(id,ct),check.names=FALSE,showProgress=FALSE)
  setnames(q,c(id,ct),c('transcript_id','estimated_count'))
  q[,transcript_id:=as.character(transcript_id)]
  q[,estimated_count:=as.numeric(estimated_count)]
  if (anyNA(q$transcript_id) || anyDuplicated(q$transcript_id) ||
      anyNA(q$estimated_count) || any(!is.finite(q$estimated_count)) ||
      any(q$estimated_count<0)) stop('Invalid Oarfish quantification: ',path)
  q
}
find_quant <- function() {
  ff <- list.files(PREMERGE_DIR,recursive=TRUE,full.names=TRUE,
                   pattern='\\.CWOWSamplesCohort_merged\\.quant$')
  ff <- ff[file.exists(ff) & !dir.exists(ff)]
  # Also accept .quant directories with a tabular quant file.
  dd <- list.dirs(PREMERGE_DIR,recursive=TRUE,full.names=TRUE)
  dd <- dd[grepl('\\.CWOWSamplesCohort_merged\\.quant$',dd)]
  for (d in dd) {
    candidates <- file.path(d,c('quant.sf','quant.tsv','abundance.tsv'))
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates)) ff <- c(ff,candidates[1])
  }
  ff <- unique(ff)
  if (!length(ff)) stop('No completed IsoMatch Oarfish .quant files found in ',PREMERGE_DIR)
  nm <- basename(ff)
  nested <- nm %in% c('quant.sf','quant.tsv','abundance.tsv')
  nm[nested] <- basename(dirname(ff[nested]))
  sample <- sub('\\.CWOWSamplesCohort_merged\\.quant$','',nm)
  z <- data.table(sample_id=sample,path=ff,size=file.info(ff)$size,
                  modified=file.info(ff)$mtime)
  z <- z[!is.na(size) & size>0]
  setorder(z,sample_id,-modified)
  if (anyDuplicated(z$sample_id)) {
    log_msg('WARNING duplicate quant files: newest selected for each sample')
    z <- unique(z,by='sample_id')
  }
  z
}
check_design <- function(md,form) {
  mm <- model.matrix(form,data=md)
  if (qr(mm)$rank<ncol(mm)) stop('Rank-deficient model: ',qr(mm)$rank,'/',ncol(mm),
                                 '. Partial cohort may confound run_pair with disease/sex.')
  if (nrow(mm)<=ncol(mm)) stop('Insufficient residual degrees of freedom')
  log_msg('Model: ',nrow(mm),' samples; ',ncol(mm),' coefficients; rank=',qr(mm)$rank)
}
# Save fitted object BEFORE contrasts so an interrupted extraction is recoverable.
fit_or_load <- function(count_mat,md,formula,path,label) {
  if (file.exists(path)) {
    log_msg('Reloading checkpoint for ',label,': ',path)
    return(readRDS(path))
  }
  dds <- timed(paste0('Construct ',label,' DESeqDataSet'),
               DESeqDataSetFromMatrix(countData=count_mat,colData=md,design=formula))
  dds <- timed(paste0('Fit ',label,' DESeq2'),DESeq(dds,quiet=TRUE,parallel=FALSE))
  timed(paste0('Save ',label,' fitted checkpoint'),saveRDS(dds,path,compress=FALSE))
  dds
}
extract_contrast <- function(dds,contrast,nm,ids,genes,md,analysis) {
  od <- file.path(OUT,nm)
  dir.create(od,recursive=TRUE,showWarnings=FALSE)
  final <- file.path(od,'DTE_results.tsv.gz')
  if (file.exists(final) && file.exists(file.path(od,'DTE_FDR05.tsv.gz')) &&
      file.exists(file.path(od,'contrast_summary.tsv'))) {
    log_msg('SKIP completed contrast ',nm)
    return(fread(file.path(od,'contrast_summary.tsv')))
  }
  chunks <- file.path(od,'contrast_chunks')
  dir.create(chunks,showWarnings=FALSE)
  starts <- seq.int(1L,length(ids),by=RESULTS_CHUNK_SIZE)
  paths <- character(length(starts))
  for (k in seq_along(starts)) {
    lo <- starts[k]; hi <- min(lo+RESULTS_CHUNK_SIZE-1L,length(ids))
    paths[k] <- file.path(chunks,sprintf('chunk_%05d.rds',k))
    if (file.exists(paths[k])) next
    log_msg('START ',nm,' chunk ',k,'/',length(starts),' rows ',lo,'-',hi)
    r <- results(dds[lo:hi,,drop=FALSE],contrast=contrast,
                 independentFiltering=FALSE,cooksCutoff=FALSE,alpha=FDR)
    x <- as.data.table(as.data.frame(r),keep.rownames='transcript_id')
    if (!identical(x$transcript_id,ids[lo:hi])) stop('Transcript ID mismatch in ',nm,' chunk ',k)
    x[,gene_id:=genes[lo:hi]]
    saveRDS(x,paths[k],compress=FALSE)
    rm(x,r); invisible(gc())
    log_msg('END ',nm,' chunk ',k)
  }
  ans <- timed(paste0('Assemble ',nm),rbindlist(lapply(paths,readRDS)))
  if (!identical(ans$transcript_id,ids)) stop('Chunk assembly transcript ID mismatch')
  # One BH adjustment per complete contrast (never per chunk).
  ans[,padj:=p.adjust(pvalue,method='BH')]
  setcolorder(ans,c('transcript_id','gene_id',setdiff(names(ans),c('transcript_id','gene_id'))))
  setorder(ans,padj)
  sig <- sum(!is.na(ans$padj) & ans$padj<FDR)
  timed(paste0('Write ',nm),{
    fwrite(ans,final,sep='\t')
    fwrite(ans[!is.na(padj) & padj<FDR],file.path(od,'DTE_FDR05.tsv.gz'),sep='\t')
    fwrite(as.data.table(md),file.path(od,'model_metadata_used.tsv'),sep='\t')
  })
  s <- data.table(analysis=analysis,comparison=nm,n_samples=nrow(md),
                  n_transcripts_tested=nrow(ans),n_significant=sig,
                  n_padj_nonNA=sum(!is.na(ans$padj)),
                  min_estimated_count=MIN_EST_COUNTS,min_samples=MIN_SAMPLES)
  fwrite(s,file.path(od,'contrast_summary.tsv'),sep='\t')
  log_msg('COMPLETED ',nm,': ',sig,' FDR<',FDR,' out of ',nrow(ans))
  rm(ans); invisible(gc())
  s
}
run_all <- function() {
  log_msg('Oarfish IsoMatch merged, moderate confidence DTE')
  log_msg('Output: ',OUT)
  log_msg('Filter: >=',MIN_EST_COUNTS,' estimated counts in >=',MIN_SAMPLES,' completed samples')
  log_msg('NOTE: Oarfish fractional estimates rounded for DESeq2 after filtering')
  map <- timed('Read moderate transcript whitelist',fread(MAP_FILE,select=c('transcript_id','gene_id')))
  if (anyDuplicated(map$transcript_id) || anyNA(map$gene_id)) stop('Invalid moderate map')
  meta <- timed('Read metadata',fread(METADATA_FILE,check.names=FALSE))
  needed <- c('sample_id','group','sex','Age','RIN_workbook','all_pacbio_runs','TDP.43')
  if (!all(needed %in% names(meta))) stop('Missing metadata: ',paste(setdiff(needed,names(meta)),collapse=', '))
  meta[,sample_id:=trimws(as.character(sample_id))]
  meta[,group:=gsub('-','_',toupper(trimws(as.character(group))))]
  meta[,sex:=toupper(trimws(as.character(sex)))]
  meta[sex=='MALE',sex:='XY'];meta[sex=='FEMALE',sex:='XX']
  meta[,Age:=as.numeric(Age)]
  meta[,RIN:=as.numeric(RIN_workbook)]
  meta[,TDP43:=as.numeric(`TDP.43`)]
  meta[,run_pair_original:=vapply(as.character(all_pacbio_runs),clean_pair,character(1))]
  keys <- sort(unique(meta$run_pair_original[!is.na(meta$run_pair_original)]))
  run_map <- data.table(run_pair_original=keys,run_pair=sprintf('RunPair_%02d',seq_along(keys)))
  meta[run_map,on='run_pair_original',run_pair:=i.run_pair]
  fwrite(run_map,file.path(OUT,'run_pair_safe_label_mapping.tsv'),sep='\t')
  if (anyDuplicated(meta$sample_id)) stop('Duplicate metadata IDs')
  if (any(!meta$group %in% GROUPS) || any(!meta$sex %in% c('XX','XY'))) stop('Invalid group/sex metadata')
  manifest <- timed('Discover Oarfish quant files',find_quant())
  manifest <- manifest[sample_id %in% meta$sample_id]
  if (!nrow(manifest)) stop('No quant files match metadata sample IDs')
  setorder(manifest,sample_id)
  manifest[,group:=meta$group[match(sample_id,meta$sample_id)]]
  manifest[,group_label:=GROUP_LABELS[group]]
  manifest[,RIN_available:=!is.na(meta$RIN[match(sample_id,meta$sample_id)])]
  fwrite(manifest,file.path(OUT,'quant_manifest.tsv'),sep='\t')
  log_msg('Discovered ',nrow(manifest),' samples; groups: ',
          paste(names(table(manifest$group)),as.integer(table(manifest$group)),collapse=', '))
  # Read one sample at a time; retain ONLY moderate transcript IDs in the matrix.
  n <- nrow(map); m <- nrow(manifest)
  mat <- matrix(0,nrow=n,ncol=m,dimnames=list(map$transcript_id,manifest$sample_id))
  sample_stats <- vector('list',m)
  for (j in seq_len(m)) {
    log_msg('Read sample ',j,'/',m,': ',manifest$sample_id[j])
    q <- tryCatch(read_quant(manifest$path[j]),error=function(e) stop('Cannot read ',manifest$sample_id[j],': ',conditionMessage(e)))
    idx <- match(q$transcript_id,map$transcript_id)
    ok <- !is.na(idx)
    if (!any(ok)) stop('ZERO transcript IDs match moderate map for ',manifest$sample_id[j],
                       '; inspect Oarfish reference identifiers')
    mat[idx[ok],j] <- q$estimated_count[ok]
    sample_stats[[j]] <- data.table(sample_id=manifest$sample_id[j],
      n_quant_transcripts=nrow(q),n_moderate_matched=sum(ok),
      sum_estimated_all=sum(q$estimated_count),sum_estimated_moderate=sum(q$estimated_count[ok]))
    rm(q,idx,ok); if (j%%5L==0L) {invisible(gc());mem(paste0('After sample ',j))}
  }
  fwrite(rbindlist(sample_stats),file.path(OUT,'per_sample_quantification_summary.tsv'),sep='\t')
  # Global prevalence is based on all completed quant files, including any missing RIN.
  present <- integer(n)
  for (j in seq_len(m)) present <- present+as.integer(mat[,j]>=MIN_EST_COUNTS)
  keep <- present>=MIN_SAMPLES
  n_keep <- sum(keep)
  log_msg('Retained ',n_keep,' / ',n,' moderate transcripts using ',m,' completed samples')
  fwrite(data.table(transcript_id=map$transcript_id[keep],gene_id=map$gene_id[keep],
                    n_samples_ge_cutoff=present[keep]),
         file.path(OUT,'globally_retained_transcripts.tsv.gz'),sep='\t')
  fwrite(data.table(n_completed_samples=m,n_moderate_transcripts=n,
                    n_retained=n_keep,min_estimated_count=MIN_EST_COUNTS,
                    min_samples=MIN_SAMPLES),
         file.path(OUT,'global_filter_summary.tsv'),sep='\t')
  if (n_keep==0L) stop('No transcripts passed global filter')
  if (n_keep!=171084L) log_msg('NOTE expected 171084 from earlier sensitivity analysis, but this snapshot retains ',
                               n_keep,'. Sample availability may have changed.')
  ids <- map$transcript_id[keep]; genes <- map$gene_id[keep]
  mat <- mat[keep,,drop=FALSE]; rm(map,present,keep); invisible(gc())
  if (any(mat>.Machine$integer.max)) stop('Oarfish estimated count exceeds R integer limit')
  # Explicitly record the rounding approximation; no TPM/length-scaled values used.
  mat <- round(mat); storage.mode(mat) <- 'integer'
  meta <- meta[match(colnames(mat),sample_id)]
  model_ix <- which(!is.na(meta$Age) & !is.na(meta$RIN) & !is.na(meta$run_pair))
  md <- as.data.frame(meta[model_ix]); rownames(md) <- md$sample_id
  md$group <- factor(md$group,levels=GROUPS)
  md$sex <- droplevels(factor(md$sex))
  md$run_pair <- droplevels(factor(md$run_pair))
  md$Age_centered <- (md$Age-mean(md$Age))/10
  md$RIN <- md$RIN-mean(md$RIN)
  group_n <- table(md$group)
  fwrite(data.table(group=names(group_n),label=GROUP_LABELS[names(group_n)],n=as.integer(group_n)),
         file.path(OUT,'disease_model_group_sizes.tsv'),sep='\t')
  if (any(group_n<MIN_GROUP_N)) stop('Partial cohort lacks >=',MIN_GROUP_N,
                                     ' samples in every disease group; DTE deferred')
  form <- ~0+group+Age_centered+RIN+run_pair+sex
  check_design(md,form)
  fwrite(as.data.table(md),file.path(OUT,'disease_model_metadata.tsv'),sep='\t')
  fwrite(data.table(sample_id=meta$sample_id[!seq_len(nrow(meta)) %in% model_ix]),
         file.path(OUT,'excluded_missing_model_covariates.tsv'),sep='\t')
  disease_counts <- mat[,model_ix,drop=FALSE]
  dds <- fit_or_load(disease_counts,md,form,file.path(OUT,'fitted_disease_dds.rds'),'disease')
  rm(disease_counts); invisible(gc());mem('After disease fit')
  fwrite(data.table(sample_id=colnames(dds),size_factor=sizeFactors(dds)),
         file.path(OUT,'disease_size_factors.tsv'),sep='\t')
  log_msg('Disease coefficient names: ',paste(resultsNames(dds),collapse=', '))
  summaries <- list()
  for (nm in names(comparisons)) {
    pair <- comparisons[[nm]]
    co <- resultsNames(dds)
    g1 <- paste0('group',pair[1]);g2 <- paste0('group',pair[2])
    if (!all(c(g1,g2) %in% co)) stop('Missing coefficients for ',nm)
    contrast <- setNames(numeric(length(co)),co)
    contrast[g1] <- 1;contrast[g2] <- -1
    pair_md <- md[md$group %in% pair,,drop=FALSE]
    summaries[[nm]] <- extract_contrast(dds,unname(contrast),nm,ids,genes,pair_md,'Disease_single_model')
  }
  rm(dds); invisible(gc());mem('After disease contrasts')
  if (RUN_TDP43) {
    ti <- which(!is.na(meta$TDP43) & !is.na(meta$Age) & !is.na(meta$RIN) & !is.na(meta$run_pair))
    tmd <- as.data.frame(meta[ti]);rownames(tmd) <- tmd$sample_id
    tmd$group <- droplevels(factor(tmd$group,levels=GROUPS))
    tmd$sex <- droplevels(factor(tmd$sex))
    tmd$run_pair <- droplevels(factor(tmd$run_pair))
    tmd$TDP43_status <- factor(ifelse(tmd$TDP43==1,'Positive','Negative'),
                              levels=c('Negative','Positive'))
    tmd$Age_centered <- (tmd$Age-mean(tmd$Age))/10
    tmd$RIN <- tmd$RIN-mean(tmd$RIN)
    if (any(table(tmd$TDP43_status)<MIN_GROUP_N)) {
      log_msg('TDP-43 SKIPPED: insufficient classified samples in one status group')
    } else {
      tf <- ~0+group+Age_centered+RIN+run_pair+sex+TDP43_status
      check_design(tmd,tf)
      tdds <- fit_or_load(mat[,ti,drop=FALSE],tmd,tf,
                          file.path(OUT,'fitted_TDP43_dds.rds'),'TDP43')
      summaries[['TDP43_positive_vs_negative']] <- extract_contrast(
        tdds,c('TDP43_status','Positive','Negative'),
        'TDP43_positive_vs_negative',ids,genes,tmd,'TDP43_adjusted_for_disease')
      fwrite(data.table(sample_id=colnames(tdds),size_factor=sizeFactors(tdds)),
             file.path(OUT,'TDP43_size_factors.tsv'),sep='\t')
      rm(tdds);invisible(gc())
    }
  }
  fwrite(rbindlist(summaries,use.names=TRUE,fill=TRUE),
         file.path(OUT,'DTE_all_comparisons_summary.tsv'),sep='\t')
  writeLines(c('Oarfish IsoMatch merged reference; moderate-confidence transcript whitelist',
               paste('Snapshot samples:',m),
               paste('Filter: estimated counts >=',MIN_EST_COUNTS,'in >=',MIN_SAMPLES,'samples'),
               paste('Retained transcripts:',n_keep),
               'Oarfish estimated counts rounded to nearest integer for DESeq2.',
               'Disease model: ~0+group+Age_centered+RIN+run_pair+sex',
               'TDP-43: separate model with additional TDP43_status term.',
               'BH adjusted per contrast; DESeq2 independent filtering disabled.',
               'Results are preliminary until all samples have completed.'),
             file.path(OUT,'ANALYSIS_README.txt'))
  log_msg('ALL COMPLETED: ',OUT)
}
tryCatch(run_all(),error=function(e){
  log_msg('ERROR: ',conditionMessage(e))
  stop(e)
},finally={log_msg('Log saved: ',LOG_FILE);close(LOG_CON)})
