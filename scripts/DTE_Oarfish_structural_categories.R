# CWOW | Structural categories of Oarfish IsoMatch-merged DTE transcripts
# Run with source() in Posit Workbench. Uses the EXACT retained transcript IDs
# and sample manifest from a completed Oarfish DTE snapshot (no re-filtering).
# Structural labels are inherited from the existing post-merge tracking/SQANTI3
# audit; they are NOT a de novo structural reclassification.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# -------------------------- EDIT THESE PATHS -------------------------------
DTE_ROOT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/DTE_Oarfish_IsoMatchMerge'
# NULL = automatically select latest snapshot containing the retained IDs and manifest.
# Set to a specific snapshot_* directory to reproduce the same poster analysis.
DTE_SNAPSHOT <- file.path(DTE_ROOT, 'snapshot_20261012_174918')
# Output from 03_PostMerge_Moderate_Structural_Category_QC.R
STRUCTURAL_ANNOTATION <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/postmerge_moderate_structural_QC/retained_transcript_structural_categories.tsv'
# Fallback location if your post-merge QC output was saved under results instead:
STRUCTURAL_ANNOTATION_ALT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs/05_DTE_DTU/postmerge_moderate_structural_QC/retained_transcript_structural_categories.tsv'
DPI <- 300
# ---------------------------------------------------------------------------

if (is.null(DTE_SNAPSHOT)) {
  ss <- list.dirs(DTE_ROOT, recursive=FALSE, full.names=TRUE)
  ss <- ss[file.exists(file.path(ss,'globally_retained_transcripts.tsv.gz')) &
             file.exists(file.path(ss,'quant_manifest.tsv'))]
  if (!length(ss)) stop('No completed DTE snapshots found under: ',DTE_ROOT)
  DTE_SNAPSHOT <- sort(ss, decreasing=TRUE)[1]
}
if (!file.exists(STRUCTURAL_ANNOTATION) && file.exists(STRUCTURAL_ANNOTATION_ALT))
  STRUCTURAL_ANNOTATION <- STRUCTURAL_ANNOTATION_ALT
for (p in c(STRUCTURAL_ANNOTATION,
            file.path(DTE_SNAPSHOT,'globally_retained_transcripts.tsv.gz'),
            file.path(DTE_SNAPSHOT,'quant_manifest.tsv')))
  if (!file.exists(p)) stop('Missing input: ',p,
                            '\nRun 03_PostMerge_Moderate_Structural_Category_QC.R first if structural annotations are absent.')
OUT <- file.path(DTE_SNAPSHOT,'Structural_Category_QC')
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
logfile <- file.path(OUT,'structural_category_run.log')
log_msg <- function(...) {
  s <- paste0(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),' | ',paste0(...,collapse=''))
  message(s); cat(s,'\n',file=logfile,append=TRUE)
}
save_table <- function(x,name) fwrite(x,file.path(OUT,name),sep='\t',na='NA')
save_plot <- function(p,stem,w=11,h=7) {
  ggsave(file.path(OUT,paste0(stem,'.png')),p,width=w,height=h,dpi=DPI,bg='white',limitsize=FALSE)
  ggsave(file.path(OUT,paste0(stem,'.pdf')),p,width=w,height=h,bg='white',limitsize=FALSE)
  log_msg('Saved ',stem)
}
CATS <- c('FSM','ISM','NIC','NNC','Fusion','Genic','Genic intron',
          'Antisense','Intergenic','More junctions','Other','Mixed','Unresolved')
COLORS <- c('FSM'='#4F83B6','ISM'='#E99545','NIC'='#F5C342','NNC'='#9D9FA1',
            'Fusion'='#75B1D0','Genic'='#79AA5B','Genic intron'='#B99A73',
            'Antisense'='#29496C','Intergenic'='#963B20',
            'More junctions'='#C65C91','Other'='#B8B8B8',
            'Mixed'='#7651A8','Unresolved'='#DADADA')
GROUPS <- c('CONTROL','AD','LBD_S','LBD_AS','LBD_ATS')
LABELS <- c(CONTROL='Control',AD='AD(AT)',LBD_S='LBD(S)',
            LBD_AS='LBD(AS)',LBD_ATS='LBD(ATS)')
THEME <- theme_classic(base_size=13) +
  theme(plot.title=element_text(face='bold',size=17),
        plot.subtitle=element_text(size=11),
        axis.text.x=element_text(color='black'),legend.position='right')

log_msg('DTE snapshot: ',DTE_SNAPSHOT)
keep <- fread(file.path(DTE_SNAPSHOT,'globally_retained_transcripts.tsv.gz'),
              select='transcript_id')
if (!nrow(keep) || anyNA(keep$transcript_id) || anyDuplicated(keep$transcript_id))
  stop('Invalid globally retained transcript list.')
log_msg('Exactly retained DTE transcripts: ',format(nrow(keep),big.mark=','))
if (nrow(keep)!=171084L) log_msg('NOTE: count differs from 171,084; plotting exact snapshot list.')
anno <- fread(STRUCTURAL_ANNOTATION)
if (!all(c('transcript_id','category') %in% names(anno)))
  stop('Annotation must contain transcript_id and category columns.')
if (anyDuplicated(anno$transcript_id)) stop('Duplicate transcript IDs in structural annotation.')
anno <- anno[keep,on='transcript_id']  # retain EXACT DTE transcript universe
anno[is.na(category) | !nzchar(category),category:='Unresolved']
if (any(!anno$category %in% CATS))
  stop('Unexpected categories: ',paste(unique(anno$category[!anno$category %in% CATS]),collapse=', '))
anno[,category:=factor(category,levels=CATS)]
missing <- sum(as.character(anno$category)=='Unresolved')
log_msg('Unresolved transcript assignments: ',missing)
if (missing/nrow(anno)>0.05)
  warning('More than 5% unresolved structural labels; inspect post-merge source mapping QC.')
save_table(anno,'DTE_retained_transcript_structural_categories.tsv')

# Plot 1: composition of UNIQUE transcript models retained for DTE.
cohort <- anno[,.(n_unique_transcripts=.N),by=category]
cohort[,percent:=100*n_unique_transcripts/sum(n_unique_transcripts)]
save_table(cohort,'DTE_unique_category_composition.tsv')
p_unique <- ggplot(cohort,aes(category,n_unique_transcripts,fill=category))+
  geom_col(width=.75)+scale_fill_manual(values=COLORS,drop=FALSE)+
  scale_y_continuous(labels=comma)+
  labs(title='Structural categories retained for Oarfish DTE',
       subtitle=paste0(format(nrow(anno),big.mark=','),
                       ' unique IsoMatch transcript models; moderate-confidence, >=5 estimated counts in >=10 samples'),
       x=NULL,y='Unique retained transcript models',fill='Structural category')+
  THEME+theme(axis.text.x=element_text(angle=40,hjust=1))
save_plot(p_unique,'Figure_01_DTE_unique_structural_category_counts',12,7)

p_unique_pct <- ggplot(cohort,aes(x='DTE retained',y=n_unique_transcripts,fill=category))+
  geom_col(position='fill',width=.6)+scale_fill_manual(values=COLORS,drop=FALSE)+
  scale_y_continuous(labels=percent_format())+
  labs(title='Structural composition of Oarfish DTE transcript universe',
       subtitle='Each retained merged transcript model counted exactly once',
       x=NULL,y='Unique transcript models (%)',fill='Structural category')+
  THEME
save_plot(p_unique_pct,'Figure_02_DTE_unique_category_100pct',9,7)

# Use the DTE snapshot manifest to calculate disease-group sample observations.
# This intentionally matches the DTE global filter cohort, including samples
# excluded from the fitted model for missing RIN; see model-only plot below.
manifest <- fread(file.path(DTE_SNAPSHOT,'quant_manifest.tsv'))
if (!all(c('sample_id','path','group') %in% names(manifest)))
  stop('Snapshot manifest must contain sample_id, path, group.')
if (anyDuplicated(manifest$sample_id)) stop('Duplicate sample IDs in manifest.')
manifest[,group:=gsub('-','_',toupper(as.character(group)))]
if (any(!manifest$group %in% GROUPS)) stop('Unrecognized group in manifest.')
if (any(!file.exists(manifest$path)))
  stop('Oarfish quantification files no longer exist: ',paste(manifest$sample_id[!file.exists(manifest$path)],collapse=', '))
model_file <- file.path(DTE_SNAPSHOT,'disease_model_metadata.tsv')
model_ids <- if (file.exists(model_file)) fread(model_file,select='sample_id')$sample_id else character()
manifest[,in_disease_model:=sample_id %in% model_ids]
save_table(manifest,'structural_QC_sample_manifest.tsv')
log_msg('Quantified samples=',nrow(manifest),'; fitted disease-model samples=',sum(manifest$in_disease_model))

read_quant <- function(path) {
  h <- names(fread(path,nrows=0,check.names=FALSE,showProgress=FALSE))
  id <- intersect(c('name','Name','target_id','target','transcript_id','transcript','tname'),h)[1]
  ct <- intersect(c('num_reads','NumReads','numReads','est_counts','est_count',
                    'estimated_counts','estimated_count','count','counts','read_count'),h)[1]
  if (is.na(id) || is.na(ct)) stop('Unexpected Oarfish quant columns in ',path,': ',paste(h,collapse=', '))
  z <- fread(path,select=c(id,ct),check.names=FALSE,showProgress=FALSE)
  setnames(z,c(id,ct),c('transcript_id','estimated_count'))
  z[,estimated_count:=as.numeric(estimated_count)]
  if (anyDuplicated(z$transcript_id) || anyNA(z$estimated_count) ||
      any(!is.finite(z$estimated_count)) || any(z$estimated_count<0))
    stop('Invalid Oarfish estimates: ',path)
  z
}
# Match the original filtering criterion exactly: >=5 estimated counts in >=10
# samples, using unrounded estimates. For per-sample category plots, "detected"
# means estimated count >0, NOT a reapplication of the >=5/10 filter.
parts <- vector('list',nrow(manifest))
ids <- anno$transcript_id
cats <- as.character(anno$category)
for (j in seq_len(nrow(manifest))) {
  log_msg('Processing Oarfish ',j,'/',nrow(manifest),' ',manifest$sample_id[j])
  q <- read_quant(manifest$path[j])
  idx <- match(ids,q$transcript_id)
  if (anyNA(idx)) {
    # Quantifications usually list the entire reference. Missing IDs treated as zero,
    # but recorded for transparency.
    log_msg('WARNING ',sum(is.na(idx)),' retained IDs absent from ',manifest$sample_id[j])
  }
  x <- numeric(length(ids)); present <- which(!is.na(idx))
  x[present] <- q$estimated_count[idx[present]]
  d <- data.table(category=factor(cats,levels=CATS),estimate=x)
  d <- d[,.(detected_transcripts=sum(estimate>0),
           transcripts_ge5=sum(estimate>=5),
           estimated_count_sum=sum(estimate)),by=category]
  d[,`:=`(sample_id=manifest$sample_id[j],group=manifest$group[j],
          in_disease_model=manifest$in_disease_model[j])]
  parts[[j]] <- d
  rm(q,d,x,idx)
  if (j%%10L==0L) invisible(gc())
}
sample_cat <- rbindlist(parts,use.names=TRUE)
sample_cat[,group_label:=factor(LABELS[group],levels=unname(LABELS[GROUPS]))]
save_table(sample_cat,'DTE_sample_category_observations.tsv')

make_group_plots <- function(dat,suffix,subsuffix) {
  agg <- dat[,.(detected_transcripts=sum(detected_transcripts),
                 transcripts_ge5=sum(transcripts_ge5),
                 estimated_count_sum=sum(estimated_count_sum)),
             by=.(group_label,category)]
  for (v in c('detected_transcripts','transcripts_ge5','estimated_count_sum')) {
    agg[,paste0('pct_',v):=if (sum(get(v))>0) 100*get(v)/sum(get(v)) else NA_real_,by=group_label]
  }
  save_table(agg,paste0('DTE_group_category_composition_',suffix,'.tsv'))
  metrics <- list(
    detected_transcripts=list(stem='detected_100pct',
      subtitle='Each retained ISOMT counted once per sample with estimated count > 0',
      ylabel='Detected transcript observations (%)'),
    transcripts_ge5=list(stem='ge5_100pct',
      subtitle='Each retained ISOMT counted once per sample with estimated count >= 5',
      ylabel='Transcript observations with >=5 counts (%)'),
    estimated_count_sum=list(stem='estimated_counts_100pct',
      subtitle='Summed Oarfish estimated counts; abundance-weighted structural composition',
      ylabel='Oarfish estimated counts (%)'))
  for (v in names(metrics)) {
    info <- metrics[[v]]
    # Plot the precomputed within-group percentages directly. This avoids
    # censoring large raw counts before position='fill' is evaluated.
    pct_col <- paste0('pct_', v)
    if (!any(is.finite(agg[[pct_col]]) & agg[[pct_col]] > 0)) {
      stop('No positive percentages for ', v, '; inspect DTE_group_category_composition_', suffix, '.tsv')
    }
    p <- ggplot(agg, aes(x=group_label, y=.data[[pct_col]], fill=category))+
      geom_col(width=.78, position='stack')+
      scale_fill_manual(values=COLORS,drop=FALSE)+
      scale_y_continuous(labels=label_percent(scale=1), breaks=seq(0,100,25),
                         expand=expansion(mult=c(0,0)))+
      coord_cartesian(ylim=c(0,100))+
      labs(title='Structural category composition by disease group',
           subtitle=paste0('Oarfish DTE retained transcripts (',subsuffix,'). ',info$subtitle),
           x=NULL,y=info$ylabel,fill='Structural category')+THEME
    save_plot(p,paste0('Figure_',suffix,'_',info$stem),12,7)
  }
}
make_group_plots(sample_cat,'03_all_filter_samples',paste0(nrow(manifest),' quantification samples'))
if (any(manifest$in_disease_model))
  make_group_plots(sample_cat[in_disease_model==TRUE],
                   '04_disease_model_samples',paste0(sum(manifest$in_disease_model),' modeled samples'))

# Model counts per disease group, for accurate poster legend.
group_n <- manifest[,.(n_quantification_samples=.N,
                        n_disease_model_samples=sum(in_disease_model)),by=group]
group_n[,group_label:=LABELS[group]]
save_table(group_n,'DTE_group_sample_counts.tsv')
writeLines(c('Structural labels are inherited from the post-merge IsoMatch tracking / rescued SQANTI3 classification.',
             'They are NOT newly calculated structural categories from the merged sequences.',
             paste('DTE snapshot:',DTE_SNAPSHOT),
             paste('Unique retained DTE transcripts:',nrow(anno)),
             paste('Samples with Oarfish quantifications:',nrow(manifest)),
             paste('Samples in disease DESeq2 model:',sum(manifest$in_disease_model)),
             'Group 100% stacked plots count sample-level transcript observations, not unique cohort isoforms.',
             'Unrounded Oarfish estimates are used in all category summaries.',
             'A global >=5 in >=10 filter is NOT the same as >=5 within each sample.',
             'Review Mixed and Unresolved category counts before biological interpretation.'),
           file.path(OUT,'README_structural_category_QC.txt'))
log_msg('DONE. Results: ',OUT)
