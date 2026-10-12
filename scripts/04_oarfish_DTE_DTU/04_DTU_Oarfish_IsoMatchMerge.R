# CWOW Oarfish IsoMatch-merged DTU (DRIMSeq)
# Source in Posit Workbench. Uses the SAME globally retained transcript IDs
# and sample snapshot as 04_DTE_Oarfish_IsoMatchMerge.R.
# DRIMSeq tests changes in relative transcript usage within genes, NOT DTE.
# All 10 disease contrasts from one adjusted model; separate TDP43 model.
# Requires data.table, DRIMSeq, BiocParallel, and the completed DTE filter snapshot.
# FIX v2: dmDSdata() requires base data.frame, not data.table.
# IMPORTANT: this script has not been executed on the user's HPC installation.

# ----------------------------- EDIT SETTINGS -----------------------------
DTE_ROOT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/DTE_Oarfish_IsoMatchMerge'
# Explicitly select the DTE snapshot with the 171084-transcript filter.
DTE_SNAPSHOT <- file.path(DTE_ROOT,'snapshot_20261012_174918')
PREMERGE_DIR <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results'
N_CPUS <- 8L
RUN_TDP43 <- TRUE
FDR <- 0.05
MIN_GROUP_N <- 3L
# If TRUE, discard genes with <2 transcripts from the DTE retained universe.
# No further transcript abundance/prevalence filter is imposed.
# New root (outside the DTE snapshot); existing checkpoints may be reused.
OUT_ROOT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/DTU_DRIMSeq'
LEGACY_OUT <- file.path(DTE_SNAPSHOT, 'DTU_DRIMSeq')
# DTU-specific thresholds applied AFTER the exact DTE 5-count/10-sample filter.
# These remove poorly informative genes and low-usage isoforms.
MIN_GENE_TOTAL <- 50L          # summed gene counts across all modeled samples
MIN_GENE_SAMPLES <- 10L        # samples with gene total >= 10
MIN_GENE_COUNT_PER_SAMPLE <- 10L
MIN_ISOFORM_PROP <- 0.01       # >= 1% of gene counts across modeled samples
MIN_ISOFORM_TOTAL <- 10L       # across modeled samples
PRECISION_SEED <- 20261012L
PRECISION_SUBSET <- 0.1
PRECISION_WORKERS <- 1L        # serial first: diagnose numerical failures

# ------------------------------------------------------------------------

for (p in c('data.table','DRIMSeq','BiocParallel','ggplot2','scales'))
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
OUT <- OUT_ROOT
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
log_msg('DRIMSeq version=',as.character(utils::packageVersion('DRIMSeq')),
        '; 8 workers registered (DRIMSeq parallel support varies by version).')
if (!'dmDSdata' %in% getNamespaceExports('DRIMSeq'))
  stop('Installed DRIMSeq does not export dmDSdata(); verify installed package/version.')
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
if (!file.exists(COUNT_CHECKPOINT) && file.exists(file.path(LEGACY_OUT,'DTU_counts.rds'))) {
  log_msg('Reuse legacy count checkpoint from ',LEGACY_OUT)
  COUNT_CHECKPOINT <- file.path(LEGACY_OUT,'DTU_counts.rds')
}
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

# DTU-specific sensitivity grid. These are ADDITIONAL filters, distinct from
# the DTE global >=5 estimated counts in >=10 samples.
# Uses integer rounded estimated counts, exactly as the fitted DRIMSeq model.
plot_filter_sensitivity <- function(wide, sample_ids, out) {
  log_msg('START DTU filter sensitivity plots')
  requireNamespace('ggplot2', quietly=TRUE) || stop('Install ggplot2')
  x <- as.matrix(wide[, ..sample_ids]); storage.mode(x) <- 'double'
  gene <- as.character(wide$gene_id)
  gene_levels <- unique(gene)
  gi <- match(gene,gene_levels)
  gene_mat <- rowsum(x,group=gi,reorder=TRUE)
  # rowsum groups sorted numeric 1..N
  gene_tot <- rowSums(gene_mat)
  n_iso <- tabulate(gi,nbins=length(gene_levels))
  grid <- CJ(min_gene_count=c(1,5,10,20), min_gene_samples=c(2,5,10,20),
             min_isoform_prop=c(0,0.005,0.01,0.02,0.05))
  iso_total <- rowSums(x)
  prop <- iso_total/pmax(gene_tot[gi],1)
  rows <- vector('list',nrow(grid))
  for (k in seq_len(nrow(grid))) {
    g_ok <- rowSums(gene_mat >= grid$min_gene_count[k]) >= grid$min_gene_samples[k]
    i_ok <- g_ok[gi] & iso_total >= MIN_ISOFORM_TOTAL & prop >= grid$min_isoform_prop[k]
    # Require at least 2 surviving isoforms per gene
    cnt <- tabulate(gi[i_ok],nbins=length(gene_levels))
    i_ok <- i_ok & cnt[gi]>=2L
    rows[[k]] <- data.table(min_gene_count=grid$min_gene_count[k],
      min_gene_samples=grid$min_gene_samples[k],
      min_isoform_prop=grid$min_isoform_prop[k],
      retained_genes=sum(cnt>=2L),retained_transcripts=sum(i_ok))
  }
  tab <- rbindlist(rows)
  fwrite(tab,file.path(out,'DTU_filter_sensitivity.tsv'),sep='\t')
  tab[,sample_threshold:=factor(min_gene_samples,levels=c(2,5,10,20))]
  for (prop0 in c(0,0.01,0.05)) {
    z <- tab[abs(min_isoform_prop-prop0)<1e-8]
    fig <- ggplot2::ggplot(z,ggplot2::aes(x=factor(min_gene_count),y=sample_threshold,fill=retained_transcripts))+
      ggplot2::geom_tile(color='white',linewidth=0.7)+
      ggplot2::geom_text(ggplot2::aes(label=format(retained_transcripts,big.mark=',')),size=3.6)+
      ggplot2::scale_fill_gradient(low='#E7F1F7',high='#27628E')+
      ggplot2::labs(title=paste0('DTU transcript retention | minimum isoform proportion ',100*prop0,'%'),
        subtitle='After DTE global filtering; genes must retain >=2 isoforms',
        x='Minimum gene counts per sample',y='Minimum samples',fill='Transcripts')+
      ggplot2::theme_minimal(base_size=12)
    ggplot2::ggsave(file.path(out,sprintf('DTU_retention_heatmap_prop_%03d.png',round(prop0*1000))),fig,width=10,height=6,dpi=300)
    ggplot2::ggsave(file.path(out,sprintf('DTU_retention_heatmap_prop_%03d.pdf',round(prop0*1000))),fig,width=10,height=6)
  }
  for (metric in c('retained_transcripts','retained_genes')) {
    fig <- ggplot2::ggplot(tab[min_isoform_prop %in% c(0,0.01,0.05)],
      ggplot2::aes(x=min_gene_samples,y=.data[[metric]],color=factor(min_gene_count),
                   linetype=factor(min_isoform_prop),group=interaction(min_gene_count,min_isoform_prop)))+
      ggplot2::geom_line(linewidth=0.8)+ggplot2::geom_point(size=2)+
      ggplot2::scale_y_continuous(labels=scales::comma)+
      ggplot2::labs(title=paste('DTU filter sensitivity:',gsub('_',' ',metric)),
        x='Minimum number of samples',y=gsub('_',' ',metric),
        color='Gene counts/sample',linetype='Min isoform proportion')+
      ggplot2::theme_minimal(base_size=12)
    ggplot2::ggsave(file.path(out,paste0('DTU_retention_curve_',metric,'.png')),fig,width=10,height=6,dpi=300)
    ggplot2::ggsave(file.path(out,paste0('DTU_retention_curve_',metric,'.pdf')),fig,width=10,height=6)
  }
  rm(x,gene_mat);invisible(gc())
  log_msg('END DTU filter sensitivity plots')
}
SENSITIVITY_DONE <- file.path(OUT,'DTU_filter_sensitivity.tsv')
if (!file.exists(SENSITIVITY_DONE)) plot_filter_sensitivity(wide,md$sample_id,OUT)

# Apply DTU-specific gene-expression and isoform-proportion filters.
# Keep this as a separate immutable checkpoint, preventing accidental reuse of
# the old unfiltered constructed dmDS from the previous script.
DTU_FILTERED_PATH <- file.path(OUT,'DTU_filtered_counts_v3.rds')
if (file.exists(DTU_FILTERED_PATH)) {
  log_msg('Load DTU-filtered count checkpoint')
  wide_filtered <- readRDS(DTU_FILTERED_PATH)
} else {
  log_msg('START DTU-specific filtering')
  model_ids <- as.character(md$sample_id)
  x <- as.matrix(wide[, ..model_ids]);storage.mode(x)<-'double'
  gi <- match(wide$gene_id,unique(wide$gene_id))
  gm <- rowsum(x,group=gi,reorder=TRUE)
  gt <- rowSums(gm);it <- rowSums(x)
  gene_ok <- gt>=MIN_GENE_TOTAL & rowSums(gm>=MIN_GENE_COUNT_PER_SAMPLE)>=MIN_GENE_SAMPLES
  iso_ok <- gene_ok[gi] & it>=MIN_ISOFORM_TOTAL & it/pmax(gt[gi],1)>=MIN_ISOFORM_PROP
  n_iso <- tabulate(gi[iso_ok],nbins=length(gt))
  keep <- iso_ok & n_iso[gi]>=2L
  wide_filtered <- wide[keep]
  log_msg('DTU filter retained ',nrow(wide_filtered),' transcripts in ',uniqueN(wide_filtered$gene_id),' genes')
  if (!nrow(wide_filtered)) stop('DTU filters retained no transcripts')
  atomic_rds(wide_filtered,DTU_FILTERED_PATH)
  fwrite(wide_filtered[,.(gene_id,feature_id)],file.path(OUT,'DTU_filtered_transcripts_v3.tsv.gz'),sep='\t')
  rm(x,gm);invisible(gc())
}
wide <- wide_filtered;rm(wide_filtered);invisible(gc())
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
  if (file.exists(path)) {
    log_msg('Load fitted ',label,' checkpoint')
    return(readRDS(path))
  }
  if (!all(meta$sample_id %in% names(wide))) stop('Wide counts missing model samples')
  # DRIMSeq's exported constructor is dmDSdata(), NOT dmDS().
  # The original script stopped before fitting because dmDS is not exported.
  # Save intermediate checkpoints to avoid repeating construction/precision.
  construct_path <- file.path(OUT,paste0('constructed_',label,'_dmDS_v3.rds'))
  precision_path <- file.path(OUT,paste0('precision_',label,'_dmDS_v3.rds'))
  if (file.exists(precision_path)) {
    log_msg('Load precision ',label,' checkpoint')
    d <- readRDS(precision_path)
  } else {
    if (file.exists(construct_path)) {
      log_msg('Load constructed ',label,' checkpoint')
      d <- readRDS(construct_path)
    } else {
      log_msg('Before construction ',label,'; R GC used MB=',
              round(sum(gc()[,2]),1))
      # dmDSdata requires class(counts) == 'data.frame' EXACTLY.
      # data.table inherits data.frame but has class c('data.table','data.frame'),
      # so passing a data.table fails the constructor's stopifnot().
      # Convert the selected columns explicitly; do not convert the full
      # checkpoint before selecting the model's samples.
      x <- as.data.frame(wide[,c('gene_id','feature_id',meta$sample_id),with=FALSE])
      stopifnot(identical(class(x), 'data.frame'))
      samples_df <- data.frame(sample_id=as.character(meta$sample_id),
                               stringsAsFactors=FALSE)
      stopifnot(identical(names(x)[-(1:2)], samples_df$sample_id))
      log_msg('dmDS input class=',paste(class(x),collapse=','),
              '; rows=',nrow(x),'; sample columns=',nrow(samples_df))
      d <- timed(paste0('Construct ',label,' dmDS'),
                 DRIMSeq::dmDSdata(counts=x, samples=samples_df))
      rm(x);invisible(gc())
      timed(paste0('Save constructed ',label,' dmDS'),
            atomic_rds(d,construct_path))
    }
    # DRIMSeq precision is estimated serially to avoid opaque worker errors.
    # Use the documented subset and seed, and fail with diagnostic information.
    set.seed(PRECISION_SEED)
    log_msg('Precision settings: workers=',PRECISION_WORKERS,
            '; subset=',PRECISION_SUBSET,'; seed=',PRECISION_SEED)
    serial_bp <- BiocParallel::SerialParam()
    d <- tryCatch(timed(paste0('Estimate ',label,' precision'),
      DRIMSeq::dmPrecision(d,design=design, BPPARAM=serial_bp)),error=function(e) {
        log_msg('PRECISION ERROR: ',conditionMessage(e))
        stop('Precision failed after DTU filtering. Inspect diagnostic counts and warnings. ',
             'Try MIN_GENE_TOTAL=100 or 250 and remove v3 checkpoints before rerun. Error: ',
             conditionMessage(e))
      })
    timed(paste0('Save precision ',label,' dmDS'),
          atomic_rds(d,precision_path))
  }
  d <- timed(paste0('Fit ',label,' DM regression'),
             DRIMSeq::dmFit(d,design=design,BPPARAM=BP))
  timed(paste0('Save ',label,' fit'),atomic_rds(d,path))
  d
}

d <- fit_model(wide,md,mm,file.path(OUT,'fitted_disease_dmDS_v3.rds'),'disease')
summary_rows <- list()
for (nm in names(comparisons)) {
  dest <- file.path(OUT,nm);dir.create(dest,showWarnings=FALSE)
  gp <- file.path(dest,'DTU_gene_results_v3.tsv.gz')
  fp <- file.path(dest,'DTU_transcript_results_v3.tsv.gz')
  sp <- file.path(dest,'DTU_summary_v3.tsv')
  if (all(file.exists(c(gp,fp,sp)))) {
    log_msg('SKIP completed ',nm);summary_rows[[nm]] <- fread(sp);next
  }
  pair <- comparisons[[nm]]
  cn <- colnames(mm); contrast <- setNames(numeric(length(cn)),cn)
  a <- paste0('group',pair[1]);b <- paste0('group',pair[2])
  if (!all(c(a,b) %in% cn)) stop('Missing disease coefficients ',nm)
  contrast[a] <- 1;contrast[b] <- -1
  test_path <- file.path(dest,'tested_dmDS_v3.rds')
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
    td <- fit_model(wide,tmd,tm,file.path(OUT,'fitted_TDP43_dmDS_v3.rds'),'TDP43')
    nm <- 'TDP43_positive_vs_negative';dest <- file.path(OUT,nm)
    dir.create(dest,showWarnings=FALSE)
    gp <- file.path(dest,'DTU_gene_results_v3.tsv.gz')
    fp <- file.path(dest,'DTU_transcript_results_v3.tsv.gz')
    sp <- file.path(dest,'DTU_summary_v3.tsv')
    if (!all(file.exists(c(gp,fp,sp)))) {
      test_path <- file.path(dest,'tested_dmDS_v3.rds')
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
fwrite(rbindlist(summary_rows,fill=TRUE),file.path(OUT,'DTU_all_comparisons_summary_v3.tsv'),sep='\t')
writeLines(c('DRIMSeq differential transcript usage (within-gene transcript proportions).',
             paste('DTE snapshot:',DTE_SNAPSHOT),
             'Starts with DTE global >=5 estimated counts in >=10 samples, then applies DTU-specific filters.',
             'DTU-specific: minimum gene total, prevalence, isoform proportion and total counts; >=2 surviving isoforms per gene.',
             'Estimated Oarfish counts rounded to integer, as in DTE.',
             'Disease: ~0+group+Age_centered+RIN+run_pair+sex; 10 pairwise contrasts.',
             'TDP43: separate adjusted model with TDP43_status; classified samples only.',
             'BH adjustment per contrast; gene-level significance primary.',
             'Partial-cohort results are preliminary. DRIMSeq requires independent validation.'),
           file.path(OUT,'DTU_ANALYSIS_README.txt'))
log_msg('ALL COMPLETE; output=',OUT)
