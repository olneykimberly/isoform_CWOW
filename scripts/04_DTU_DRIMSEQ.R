# CWOW Oarfish IsoMatch-merged DTU (DRIMSeq)
# Source in Posit Workbench. Uses the SAME globally retained transcript IDs
# and sample snapshot as 04_DTE_Oarfish_IsoMatchMerge.R.
# DRIMSeq tests changes in relative transcript usage within genes, NOT DTE.
# All 10 disease contrasts from one adjusted model; separate TDP43 model.
# Requires data.table, DRIMSeq, BiocParallel, and the completed DTE filter snapshot.
# IMPORTANT: this script has not been executed on the user's HPC installation.
# ----------------------------- EDIT SETTINGS -----------------------------
DTE_ROOT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/DTE_Oarfish_IsoMatchMerge'
# Explicitly select the DTE snapshot with the 171084-transcript filter.
DTE_SNAPSHOT <- NULL  # e.g. file.path(DTE_ROOT,'snapshot_20261012_123456')
PREMERGE_DIR <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results'
N_CPUS <- 8L
RUN_TDP43 <- TRUE
FDR <- 0.05
MIN_GROUP_N <- 3L
# If TRUE, discard genes with <2 transcripts from the DTE retained universe.
# No further transcript abundance/prevalence filter is imposed.
OUT_SUBDIR <- 'DTU_DRIMSeq'
# ------------------------------------------------------------------------

for (p in c('data.table','DRIMSeq','BiocParallel'))
  if (!requireNamespace(p,quietly=TRUE)) stop('Install required R package: ',p)
suppressPackageStartupMessages({library(data.table);library(DRIMSeq);library(BiocParallel)})
if (.Platform$OS.type != 'unix') stop('MulticoreParam requires Linux/Unix')
if (is.null(DTE_SNAPSHOT)) {
  candidates <- list.dirs(DTE_ROOT,recursive=FALSE,full.names=TRUE)
  candidates <- candidates[file.exists(file.path(candidates,'globally_retained_transcripts.tsv.gz')) &
                             file.exists(file.path(candidates,'quant_manifest.tsv')) &
                             file.exists(file.path(candidates,'disease_model_metadata.tsv'))]
  if (!length(candidates)) stop('No completed DTE snapshot with filter/manifest/model metadata found. Set DTE_SNAPSHOT.')
  DTE_SNAPSHOT <- sort(candidates,decreasing=TRUE)[1]
}
for (f in c('globally_retained_transcripts.tsv.gz','quant_manifest.tsv','disease_model_metadata.tsv'))
  if (!file.exists(file.path(DTE_SNAPSHOT,f))) stop('Missing DTE snapshot file: ',f)
OUT <- file.path(DTE_SNAPSHOT,OUT_SUBDIR)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
LOG <- file.path(OUT,'DTU_DRIMSeq.log')
log_msg <- function(...) {
  z <- paste(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),paste0(...,collapse=''))
  cat(z,'\n');cat(z,'\n',file=LOG,append=TRUE)
}
timed <- function(label,expr) {
  log_msg('START ',label);t <- proc.time();x <- force(expr)
  dt <- proc.time()-t
  log_msg('END ',label,' wall_s=',round(unname(dt['elapsed']),1),
          ' parent_user_cpu_s=',round(unname(dt['user.self']),1),
          ' parent_system_cpu_s=',round(unname(dt['sys.self']),1));x
}
atomic_rds <- function(x,path) {
  tmp <- paste0(path,'.writing');saveRDS(x,tmp,compress=FALSE)
  if (!file.rename(tmp,path)) stop('Cannot finalize checkpoint: ',path)
}
read_quant <- function(path) {
  h <- names(fread(path,nrows=0,showProgress=FALSE))
  id <- intersect(c('name','Name','target_id','target','transcript_id','transcript','tname'),h)[1]
  ct <- intersect(c('num_reads','NumReads','numReads','est_counts','est_count',
                    'estimated_counts','estimated_count','count','counts','read_count'),h)[1]
  if (is.na(id)||is.na(ct)) stop('Unrecognized Oarfish quantification columns: ',path,
                                ' columns=',paste(h,collapse=','))
  z <- fread(path,select=c(id,ct),showProgress=FALSE)
  setnames(z,c(id,ct),c('transcript_id','estimated_count'))
  if (anyDuplicated(z$transcript_id) || anyNA(z$estimated_count) ||
      any(!is.finite(z$estimated_count)) || any(z$estimated_count<0))
    stop('Invalid quantification: ',path)
  z
}
GROUPS <- c('CONTROL','AD','LBD_S','LBD_AS','LBD_ATS')
comparisons <- list(
 AD_vs_CONTROL=c('AD','CONTROL'), LBD_S_vs_CONTROL=c('LBD_S','CONTROL'),
 LBD_AS_vs_CONTROL=c('LBD_AS','CONTROL'), LBD_ATS_vs_CONTROL=c('LBD_ATS','CONTROL'),
 LBD_S_vs_AD=c('LBD_S','AD'), LBD_AS_vs_AD=c('LBD_AS','AD'),
 LBD_ATS_vs_AD=c('LBD_ATS','AD'), LBD_AS_vs_LBD_S=c('LBD_AS','LBD_S'),
 LBD_ATS_vs_LBD_S=c('LBD_ATS','LBD_S'), LBD_ATS_vs_LBD_AS=c('LBD_ATS','LBD_AS'))
BP <- MulticoreParam(workers=N_CPUS,progressbar=FALSE)
allocated <- suppressWarnings(as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','')))
if (!is.na(allocated) && allocated<N_CPUS) stop('SLURM allocation smaller than N_CPUS')
log_msg('DTE snapshot=',DTE_SNAPSHOT,'; requested CPUs=',N_CPUS)
log_msg('DRIMSeq may have version-specific parallel behavior; 8 workers are registered.')
BiocParallel::register(BP,default=TRUE)

# Load the exact globally filtered IDs and gene assignments used by DTE.
retained <- fread(file.path(DTE_SNAPSHOT,'globally_retained_transcripts.tsv.gz'))
stopifnot(all(c('transcript_id','gene_id') %in% names(retained)))
if (anyDuplicated(retained$transcript_id) || anyNA(retained$gene_id)) stop('Invalid retained transcript map')
filter_summary <- fread(file.path(DTE_SNAPSHOT,'global_filter_summary.tsv'))
if (!all(c('min_estimated_count','min_samples') %in% names(filter_summary)))
  stop('Missing filter parameters')
if (filter_summary$min_estimated_count[1]!=5 || filter_summary$min_samples[1]!=10)
  stop('DTE snapshot did not use >=5 estimated counts in >=10 samples')
log_msg('DTE global filter retained ',nrow(retained),' transcripts (expected ~171084)')
manifest <- fread(file.path(DTE_SNAPSHOT,'quant_manifest.tsv'))
md <- fread(file.path(DTE_SNAPSHOT,'disease_model_metadata.tsv'))
if (!all(c('sample_id','path') %in% names(manifest))) stop('Invalid manifest')
if (!all(c('sample_id','group','Age_centered','RIN','run_pair','sex') %in% names(md)))
  stop('Invalid disease model metadata')
if (anyDuplicated(manifest$sample_id)||anyDuplicated(md$sample_id)) stop('Duplicate samples')
if (!all(md$sample_id %in% manifest$sample_id)) stop('Model sample absent from quant manifest')
# Use ALL available samples for filtering (already performed in DTE); DTU model
# uses precisely the DTE disease-model samples, not a reselected cohort.
manifest <- manifest[match(md$sample_id,manifest$sample_id)]
if (anyNA(manifest$path)) stop('Missing quant path for model sample')

# DRIMSeq only supports usage comparisons for genes with >=2 retained isoforms.
retained <- retained[,if (.N>=2L) .SD,by=gene_id]
setorder(retained,gene_id,transcript_id)
log_msg('DTU-eligible: ',nrow(retained),' transcripts in ',uniqueN(retained$gene_id),' multi-isoform genes')
fwrite(retained,file.path(OUT,'DTU_eligible_transcripts.tsv.gz'),sep='\t')
# Store sparse-ish long count table; zeros for all sample/isoform pairs are
# explicitly represented because DRIMSeq expects complete count vectors.
COUNT_CHECKPOINT <- file.path(OUT,'DTU_counts.rds')
if (file.exists(COUNT_CHECKPOINT)) {
  log_msg('Load counts checkpoint');counts <- readRDS(COUNT_CHECKPOINT)
} else {
  ids <- retained$transcript_id
  n <- length(ids);m <- nrow(manifest)
  mat <- matrix(0L,nrow=n,ncol=m,dimnames=list(ids,manifest$sample_id))
  for (j in seq_len(m)) {
    log_msg('Read quant ',j,'/',m,': ',manifest$sample_id[j])
    q <- read_quant(manifest$path[j]);ix <- match(q$transcript_id,ids)
    ok <- !is.na(ix)
    if (!any(ok)) stop('No transcript IDs matched for ',manifest$sample_id[j])
    vals <- q$estimated_count[ok]
    if (any(vals>.Machine$integer.max)) stop('Estimated count exceeds integer limit')
    mat[ix[ok],j] <- as.integer(round(vals))
    rm(q,ix,ok,vals)
    if (j%%10L==0L) invisible(gc())
  }
  # Same explicit rounding approximation as the DTE script.
  counts <- data.table(gene_id=retained$gene_id,feature_id=ids)
  for (j in seq_len(m)) counts[[colnames(mat)[j]]] <- mat[,j]
  rm(mat);invisible(gc())
  timed('Save count checkpoint',atomic_rds(counts,COUNT_CHECKPOINT))
}
if (!identical(sort(setdiff(names(counts),c('gene_id','feature_id'))),sort(md$sample_id)))
  stop('Count checkpoint samples do not match disease model snapshot')
# DRIMSeq requires counts wide: gene_id, feature_id, sample count columns.
wide <- counts
rm(counts);invisible(gc())
md <- as.data.frame(md)
md <- md[match(manifest$sample_id,md$sample_id),,drop=FALSE]
rownames(md) <- md$sample_id
md$group <- factor(as.character(md$group),levels=GROUPS)
md$sex <- droplevels(factor(md$sex))
md$run_pair <- droplevels(factor(md$run_pair))
if (any(table(md$group)<MIN_GROUP_N)) stop('Insufficient samples in disease groups')
form <- ~0+group+Age_centered+RIN+run_pair+sex
mm <- model.matrix(form,md)
if (qr(mm)$rank<ncol(mm)) stop('Rank-deficient disease design')
log_msg('Disease model samples=',nrow(md),' coefficients=',ncol(mm))
# DRIMSeq does not need DESeq2 size factors: it models within-gene proportions.
fit_model <- function(wide,meta,design,path,label) {
  if (file.exists(path)) {log_msg('Load fitted ',label,' checkpoint');return(readRDS(path))}
  if (!all(meta$sample_id %in% names(wide))) stop('Wide counts missing model samples')
  x <- wide[,c('gene_id','feature_id',meta$sample_id),with=FALSE]
  d <- timed(paste0('Construct ',label,' dmDS'),
             DRIMSeq::dmDS(counts=as.data.frame(x),samples=data.frame(sample_id=meta$sample_id)))
  # No additional expression filtering: the DTE global filter was already applied.
  # Genes with >=2 retained transcripts were selected above.
  d <- timed(paste0('Estimate ',label,' precision'),DRIMSeq::dmPrecision(d,design=design))
  d <- timed(paste0('Fit ',label,' DM regression'),DRIMSeq::dmFit(d,design=design))
  timed(paste0('Save ',label,' fit'),atomic_rds(d,path))
  d
}

d <- fit_model(wide,md,mm,file.path(OUT,'fitted_disease_dmDS.rds'),'disease')
summary_rows <- list()
for (nm in names(comparisons)) {
  dest <- file.path(OUT,nm);dir.create(dest,showWarnings=FALSE)
  gp <- file.path(dest,'DTU_gene_results.tsv.gz')
  fp <- file.path(dest,'DTU_transcript_results.tsv.gz')
  sp <- file.path(dest,'DTU_summary.tsv')
  if (all(file.exists(c(gp,fp,sp)))) {
    log_msg('SKIP completed ',nm);summary_rows[[nm]] <- fread(sp);next
  }
  pair <- comparisons[[nm]]
  cn <- colnames(mm); contrast <- setNames(numeric(length(cn)),cn)
  a <- paste0('group',pair[1]);b <- paste0('group',pair[2])
  if (!all(c(a,b) %in% cn)) stop('Missing disease coefficients ',nm)
  contrast[a] <- 1;contrast[b] <- -1
  test_path <- file.path(dest,'tested_dmDS.rds')
  if (file.exists(test_path)) {
    dt <- readRDS(test_path)
  } else {
    dt <- timed(paste0('Test ',nm),DRIMSeq::dmTest(d,contrast=unname(contrast)))
    atomic_rds(dt,test_path)
  }
  gene <- as.data.table(DRIMSeq::results(dt,level='gene'))
  feature <- as.data.table(DRIMSeq::results(dt,level='feature'))
  # BH correction separately for each contrast and each result level.
  if (!'pvalue' %in% names(gene)||!'pvalue' %in% names(feature))
    stop('Unexpected DRIMSeq result columns: ',paste(names(gene),collapse=','))
  gene[,padj:=p.adjust(pvalue,method='BH')]
  feature[,padj:=p.adjust(pvalue,method='BH')]
  fwrite(gene,gp,sep='\t');fwrite(feature,fp,sep='\t')
  s <- data.table(comparison=nm,n_samples=nrow(md),
                  n_genes_tested=nrow(gene),n_transcripts_tested=nrow(feature),
                  n_significant_genes=sum(gene$padj<FDR,na.rm=TRUE),
                  n_significant_transcripts=sum(feature$padj<FDR,na.rm=TRUE),
                  fdr=FDR)
  fwrite(s,sp,sep='\t');summary_rows[[nm]] <- s
  log_msg('COMPLETE ',nm,' significant genes=',s$n_significant_genes)
}

# TDP43: read the same DTE snapshot manifest, but include classified samples
# with complete covariates (including those not in disease model only if they
# were in the quant manifest). The existing DTE model metadata does not
# necessarily include TDP43 status, so retrieve the original metadata file.
if (RUN_TDP43) {
  META_PATH <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/QC_technical_RIN_normalization_CCA/CWOW_n274_DTE_metadata_RIN_SMRTcell_RGID_augmented.tsv'
  if (!file.exists(META_PATH)) stop('Missing original metadata for TDP43: ',META_PATH)
  orig <- fread(META_PATH)
  if (!'TDP.43' %in% names(orig)) stop('TDP.43 column missing')
  orig[,TDP43:=suppressWarnings(as.numeric(`TDP.43`))]
  # Reuse modeled disease samples to keep covariate values identical.
  tmd <- merge(md,orig[,.(sample_id,TDP43)],by='sample_id',sort=FALSE)
  tmd <- tmd[!is.na(tmd$TDP43) & tmd$TDP43 %in% c(0,1),,drop=FALSE]
  tmd$TDP43_status <- factor(ifelse(tmd$TDP43==1,'Positive','Negative'),
                             levels=c('Negative','Positive'))
  if (nrow(tmd)<10 || any(table(tmd$TDP43_status)<MIN_GROUP_N)) {
    log_msg('TDP43 SKIPPED: insufficient classified samples')
  } else {
    tmd$group <- droplevels(tmd$group);tmd$sex <- droplevels(tmd$sex)
    tmd$run_pair <- droplevels(tmd$run_pair)
    tm <- model.matrix(~0+group+Age_centered+RIN+run_pair+sex+TDP43_status,tmd)
    if (qr(tm)$rank<ncol(tm)) stop('Rank-deficient TDP43 design')
    td <- fit_model(wide,tmd,tm,file.path(OUT,'fitted_TDP43_dmDS.rds'),'TDP43')
    nm <- 'TDP43_positive_vs_negative';dest <- file.path(OUT,nm)
    dir.create(dest,showWarnings=FALSE)
    gp <- file.path(dest,'DTU_gene_results.tsv.gz')
    fp <- file.path(dest,'DTU_transcript_results.tsv.gz')
    sp <- file.path(dest,'DTU_summary.tsv')
    if (!all(file.exists(c(gp,fp,sp)))) {
      test_path <- file.path(dest,'tested_dmDS.rds')
      tt <- if (file.exists(test_path)) readRDS(test_path) else
        timed('Test TDP43',DRIMSeq::dmTest(td,coef='TDP43_statusPositive'))
      if (!file.exists(test_path)) atomic_rds(tt,test_path)
      gene <- as.data.table(DRIMSeq::results(tt,level='gene'))
      feature <- as.data.table(DRIMSeq::results(tt,level='feature'))
      gene[,padj:=p.adjust(pvalue,method='BH')]
      feature[,padj:=p.adjust(pvalue,method='BH')]
      fwrite(gene,gp,sep='\t');fwrite(feature,fp,sep='\t')
      s <- data.table(comparison=nm,n_samples=nrow(tmd),n_genes_tested=nrow(gene),
                      n_transcripts_tested=nrow(feature),
                      n_significant_genes=sum(gene$padj<FDR,na.rm=TRUE),
                      n_significant_transcripts=sum(feature$padj<FDR,na.rm=TRUE),fdr=FDR)
      fwrite(s,sp,sep='\t')
    }
    summary_rows[[nm]] <- fread(sp)
  }
}
fwrite(rbindlist(summary_rows,fill=TRUE),file.path(OUT,'DTU_all_comparisons_summary.tsv'),sep='\t')
writeLines(c('DRIMSeq differential transcript usage (within-gene transcript proportions).',
             paste('DTE snapshot:',DTE_SNAPSHOT),
             'Exactly the DTE global >=5 estimated counts in >=10 samples transcript universe.',
             'Only genes with >=2 retained isoforms are DTU testable.',
             'Estimated Oarfish counts rounded to integer, as in DTE.',
             'Disease: ~0+group+Age_centered+RIN+run_pair+sex; 10 pairwise contrasts.',
             'TDP43: separate adjusted model with TDP43_status; classified samples only.',
             'BH adjustment per contrast; gene-level significance primary.',
             'Partial-cohort results are preliminary. DRIMSeq requires independent validation.'),
           file.path(OUT,'DTU_ANALYSIS_README.txt'))
log_msg('ALL COMPLETE; output=',OUT)
