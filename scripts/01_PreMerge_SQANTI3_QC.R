# ============================================================================
# CWOW | SQANTI3 PRE-MERGE QUALITY CONTROL | 274 samples | Posit Workbench
# ============================================================================
# PURPOSE
#   Summarize rescued SQANTI3 classification files before IsoMatch merging.
#   Each sample-level transcript model is counted as an observation, NOT as a
#   unique merged IsoMatch transcript. Excludes *_sqanti3_with_sr_* files.
#
# INPUTS
#   1. IsoMatch sample manifest (expected 274 samples)
#   2. CWOW cohort metadata (disease group and biological sex)
#   3. Per-sample rescued SQANTI3 classification files
#
# OUTPUT DIRECTORY: OUT (set below)
# TABLES (.tsv)
#   sample_metadata_274.tsv                   Sample IDs, group, sex
#   group_sex_sample_counts.tsv               Group-by-sex sample counts
#   missing_rescued_classifications.tsv       Missing rescued files (if any)
#   sample_QC.tsv                              Transcript and FL support per sample
#   sample_structural_categories.tsv          Categories per sample
#   sample_transcript_length_summary.tsv      Length statistics (if available)
#   sample_quality_metrics.tsv                Available SQANTI3 quality metrics
#   group_structural_categories.tsv           Group-level category counts/percent
#   cohort_structural_categories.tsv          Cohort-level category counts/percent
#   group_yield_summary.tsv                   Per-group yield summary
#   Figure_04_pie_category_percentages.tsv    Pie-chart category counts/percent
#   group_sex_structural_categories.tsv       Group-by-sex category counts (optional)
#   SQANTI3_filter_stage_file_inventory.tsv   Initial/filtered/rescued file presence
#
# FIGURES (.png and .pdf; same stem for both formats)
#   Figure_01_sample_category_proportions      Category proportions per sample
#   Figure_02_group_category_100pct            Category composition by group
#   Figure_03_group_category_absolute          Absolute category counts by group
#   Figure_04_cohort_category_counts           Cohort counts and FL support
#   Figure_04A_cohort_observed_pie              Transcript observation pie chart
#   Figure_04B_cohort_FL_pie                    FL support pie chart
#   Figure_05_group_sample_yield               Per-sample yield by group
#   Figure_06_group_FL_weighted_category_100pct FL-weighted categories (if available)
#   Figure_07_group_sample_category_distributions Sample category distributions
#   Figure_08_group_transcript_lengths         Transcript lengths (if available)
#   Figure_09_group_quality_indicators         Quality indicators (if available)
#   Figure_10_group_sex_category_100pct        Sex-stratified composition (optional)
#
# Source directly in RStudio/Posit Workbench. No command-line arguments.
# Figure numbers are filenames only; plot titles contain no figure numbers.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# --- User settings -----------------------------------------------------------
ROOT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/'
PREMERGE_DIR <- file.path(ROOT, 'scripts/grandcanyon/results')
MANIFEST <- file.path(ROOT, 'results/isomatch_outputs/manifest_allSamples_deduplicated.txt')
METADATA <- file.path(ROOT, 'scripts/CWOW_full_long_read_cohort_n276.csv')
OUT <- file.path(ROOT, 'results/premerge_SQANTI3_QC')
EXPECTED_SAMPLES <- 274L
DPI <- 300
# If TRUE, plot sex-by-disease facets as well.
MAKE_SEX_STRATIFIED <- TRUE
# Process 274 classification files one at a time; do not load all transcripts.
# Inspired by the supplied structural-category palette.
CATEGORY_LEVELS <- c('FSM','ISM','NIC','NNC','Fusion','Genic',
                     'Genic intron','Antisense','Intergenic','More junctions','Other')
CATEGORY_COLORS <- c(
  'FSM'='#4F83B6', 'ISM'='#E99545', 'NIC'='#F5C342', 'NNC'='#9D9FA1',
  'Fusion'='#75B1D0', 'Genic'='#79AA5B', 'Genic intron'='#B99A73',
  'Antisense'='#29496C', 'Intergenic'='#963B20',
  'More junctions'='#C65C91', 'Other'='#B8B8B8'
)
GROUP_LEVELS <- c('CONTROL','AD','LBD_S','LBD_AS','LBD_ATS')
GROUP_LABELS <- c('Control','AD(AT)','LBD(S)','LBD(AS)','LBD(ATS)')
stopifnot(file.exists(MANIFEST), file.exists(METADATA), dir.exists(PREMERGE_DIR))
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
save_table <- function(x, name) fwrite(x, file.path(OUT, name), sep='\t', na='NA')
save_plot <- function(p, stem, width=12, height=7) {
  ggsave(file.path(OUT,paste0(stem,'.png')),p,width=width,height=height,dpi=DPI,bg='white',limitsize=FALSE)
  ggsave(file.path(OUT,paste0(stem,'.pdf')),p,width=width,height=height,bg='white',limitsize=FALSE)
  message('Saved: ',stem)
}
plot_theme <- theme_classic(base_size=12) + theme(
  plot.title=element_text(face='bold',size=15),
  plot.subtitle=element_text(size=10),
  axis.text.x=element_text(color='black'),
  legend.position='right',
  strip.background=element_rect(fill='#F1F1F1',colour=NA),
  strip.text=element_text(face='bold'))
scale_category <- function() scale_fill_manual(values=CATEGORY_COLORS, limits=CATEGORY_LEVELS,
                                                 drop=TRUE, name='Structural category')

# --- 1. Cohort manifest and disease metadata --------------------------------
# The IsoMatch manifest may be tab-separated sample ID + GTF, or GTF paths.
man <- fread(MANIFEST, header=FALSE, sep='\t', fill=TRUE, showProgress=FALSE)
manifest_ids <- trimws(as.character(man[[1L]]))
manifest_ids <- sub('\\.sorted\\.gtf(\\.gz)?$', '', basename(manifest_ids))
manifest_ids <- unique(manifest_ids[nzchar(manifest_ids)])
if(length(manifest_ids)!=EXPECTED_SAMPLES) stop('Manifest has ',length(manifest_ids),
  ' unique IDs; expected ',EXPECTED_SAMPLES,'. Check manifest formatting.')
meta <- fread(METADATA, showProgress=FALSE)
need <- c('Study-Subject-ID','group','Gender')
if(!all(need %in% names(meta))) stop('Metadata missing columns: ',
                                      paste(setdiff(need,names(meta)),collapse=', '))
meta[,sample_id:=trimws(as.character(`Study-Subject-ID`))]
meta[,sample_id:=ifelse(grepl('^CWOW_',sample_id),sample_id,paste0('CWOW_',sample_id))]
meta[,group:=toupper(gsub('[- ]','_',trimws(as.character(group))))]
meta[,sex:=toupper(trimws(as.character(Gender)))]
meta[sex %in% c('F','FEMALE','XX'),sex:='XX']
meta[sex %in% c('M','MALE','XY'),sex:='XY']
meta <- meta[sample_id %in% manifest_ids]
if(anyDuplicated(meta$sample_id)) stop('Duplicate metadata sample IDs after matching manifest.')
if(!setequal(meta$sample_id,manifest_ids)) stop('Missing metadata for: ',
  paste(setdiff(manifest_ids,meta$sample_id),collapse=', '))
if(any(!meta$group %in% GROUP_LEVELS)) stop('Unexpected groups: ',
  paste(unique(meta$group[!meta$group %in% GROUP_LEVELS]),collapse=', '))
meta[,group:=factor(group,levels=GROUP_LEVELS,labels=GROUP_LABELS)]
meta[,sex:=factor(sex,levels=c('XX','XY'))]
meta <- meta[,.(sample_id,group,sex)]
save_table(meta,'sample_metadata_274.tsv')
save_table(meta[,.(n_samples=.N),by=.(group,sex)],'group_sex_sample_counts.tsv')

# --- 2. Locate rescued classifications --------------------------------------
message('Finding rescued SQANTI3 classifications...')
paths <- list.files(PREMERGE_DIR,pattern='_rescued_classification\\.txt(\\.gz)?$',
                    recursive=TRUE,full.names=TRUE)
paths <- paths[!grepl('sqanti3_with_sr|with_sr',paths,ignore.case=TRUE)]
ids <- sub('_rescued_classification\\.txt(\\.gz)?$','',basename(paths))
index <- data.table(sample_id=ids,path=paths)[sample_id %in% manifest_ids]
if(anyDuplicated(index$sample_id)) stop('Multiple rescued classifications for: ',
  paste(unique(index$sample_id[duplicated(index$sample_id)]),collapse=', '))
missing <- setdiff(manifest_ids,index$sample_id)
save_table(data.table(sample_id=missing),'missing_rescued_classifications.tsv')
if(length(missing)) stop('Missing rescued classifications for ',length(missing),
                        ' samples. See missing_rescued_classifications.tsv')
setorder(index,sample_id)
first_field <- function(header, candidates) {
  hit <- candidates[candidates %in% header]
  if(length(hit)) hit[1] else NA_character_
}
normalize_class <- function(x) {
  x <- tolower(trimws(as.character(x)))
  fcase(
    x %in% c('full-splice_match','full-splice match','fsm'), 'FSM',
    x %in% c('incomplete-splice_match','incomplete-splice match','ism'), 'ISM',
    x %in% c('novel_in_catalog','novel in catalog','nic'), 'NIC',
    x %in% c('novel_not_in_catalog','novel not in catalog','nnc'), 'NNC',
    grepl('fusion',x), 'Fusion',
    x %in% c('genic','genic_genomic','genic-genomic'), 'Genic',
    x %in% c('genic_intron','genic-intron'), 'Genic intron',
    grepl('antisense',x), 'Antisense',
    grepl('intergenic',x), 'Intergenic',
    grepl('more.*junction',x), 'More junctions',
    default='Other'
  )
}
truthy <- function(x) tolower(trimws(as.character(x))) %in% c('true','yes','y','1','within','supported')
sample_parts <- vector('list',nrow(index))
category_parts <- vector('list',nrow(index))
length_parts <- vector('list',nrow(index))
quality_parts <- vector('list',nrow(index))

# --- 3. Process one sample at a time -----------------------------------------
for(i in seq_len(nrow(index))) {
  path <- index$path[i]; id <- index$sample_id[i]
  header <- names(fread(path,nrows=0,showProgress=FALSE))
  if(!all(c('isoform','structural_category') %in% header))
    stop('SQANTI3 required columns missing in ',path)
  fl <- first_field(header,c('FL','FL.BioSample_1',grep('^FL\\.',header,value=TRUE)))
  len <- first_field(header,c('length','transcript_length'))
  junction <- first_field(header,c('all_canonical','canonical','non_canonical'))
  cage <- first_field(header,c('within_CAGE_peak','within_CAGE','CAGE_peak'))
  polya <- first_field(header,c('polyA_motif','polyA_motif_found','polyA_dist','within_polyA_site'))
  selected <- unique(na.omit(c('isoform','structural_category',fl,len,junction,cage,polya)))
  d <- fread(path,select=selected,showProgress=FALSE)
  if(anyDuplicated(d$isoform)) warning('Duplicated isoform IDs: ',id)
  d[,category:=normalize_class(structural_category)]
  d[,FL_support:=if(!is.na(fl)) suppressWarnings(as.numeric(get(fl))) else NA_real_]
  if(any(d$FL_support < 0,na.rm=TRUE)) stop('Negative FL support: ',id)
  sample_parts[[i]] <- data.table(sample_id=id,transcripts=nrow(d),
    total_FL_support=if(!is.na(fl)) sum(d$FL_support,na.rm=TRUE) else NA_real_,
    median_FL_per_transcript=if(!is.na(fl)) median(d$FL_support,na.rm=TRUE) else NA_real_)
  category_parts[[i]] <- d[,.(transcripts=.N,
    FL_support=if(!is.na(fl)) sum(FL_support,na.rm=TRUE) else NA_real_),by=category][,sample_id:=id][]
  if(!is.na(len)) {
    d[,tx_length:=suppressWarnings(as.numeric(get(len)))]
    length_parts[[i]] <- d[is.finite(tx_length) & tx_length>0,
      .(n=.N,median_length=median(tx_length),q25=as.numeric(quantile(tx_length,.25)),
        q75=as.numeric(quantile(tx_length,.75))),by=category][,sample_id:=id][]
  }
  q <- list()
  if(!is.na(junction)) {
    v <- as.character(d[[junction]])
    if(junction=='non_canonical') {
      numeric_v <- suppressWarnings(as.numeric(v))
      q[['canonical_splice']] <- if(all(is.na(numeric_v))) mean(!truthy(v),na.rm=TRUE)
                                 else mean(numeric_v==0,na.rm=TRUE)
    } else q[['canonical_splice']] <- mean(truthy(v),na.rm=TRUE)
  }
  if(!is.na(cage)) q[['CAGE_peak_support']] <- mean(truthy(d[[cage]]),na.rm=TRUE)
  if(!is.na(polya)) {
    v <- as.character(d[[polya]])
    if(polya=='polyA_dist') {
      nv <- suppressWarnings(as.numeric(v))
      q[['polyA_site_within_50nt']] <- mean(is.finite(nv)&abs(nv)<=50,na.rm=TRUE)
    } else q[['polyA_evidence']] <- mean(!is.na(v)&nzchar(v)&
      !tolower(v) %in% c('na','none','false','0'),na.rm=TRUE)
  }
  if(length(q)) quality_parts[[i]] <- data.table(sample_id=id,metric=names(q),fraction=unlist(q))
  if(i %% 20L == 0L || i==nrow(index)) message('Processed ',i,'/',nrow(index),' samples')
}
samples <- merge(rbindlist(sample_parts,fill=TRUE),meta,by='sample_id',sort=FALSE)
cats <- merge(rbindlist(category_parts,fill=TRUE),meta,by='sample_id',sort=FALSE)
lengths <- rbindlist(length_parts,fill=TRUE)
if(nrow(lengths)) lengths <- merge(lengths,meta,by='sample_id',sort=FALSE)
quality <- rbindlist(quality_parts,fill=TRUE)
if(nrow(quality)) quality <- merge(quality,meta,by='sample_id',sort=FALSE)
cats[,category:=factor(category,levels=CATEGORY_LEVELS)]
if(nrow(lengths)) lengths[,category:=factor(category,levels=CATEGORY_LEVELS)]
save_table(samples,'sample_QC.tsv')
save_table(cats,'sample_structural_categories.tsv')
if(nrow(lengths)) save_table(lengths,'sample_transcript_length_summary.tsv')
if(nrow(quality)) save_table(quality,'sample_quality_metrics.tsv')

# --- 4. Summaries at sample, group and cohort levels ------------------------
group_cat <- cats[,.(transcripts=sum(transcripts),
  FL_support=if(all(is.na(FL_support))) NA_real_ else sum(FL_support,na.rm=TRUE)),
  by=.(group,category)]
group_cat[,pct_transcripts:=100*transcripts/sum(transcripts),by=group]
group_cat[,pct_FL:=if(all(is.na(FL_support))) NA_real_ else 100*FL_support/sum(FL_support,na.rm=TRUE),by=group]
cohort_cat <- cats[,.(transcripts=sum(transcripts),
  FL_support=if(all(is.na(FL_support))) NA_real_ else sum(FL_support,na.rm=TRUE)),by=category]
cohort_cat[,pct_transcripts:=100*transcripts/sum(transcripts)]
cohort_cat[,pct_FL:=if(all(is.na(FL_support))) NA_real_ else 100*FL_support/sum(FL_support,na.rm=TRUE)]
group_yield <- samples[,.(n_samples=.N,median_transcripts=median(transcripts),
  mean_transcripts=mean(transcripts),median_FL=median(total_FL_support,na.rm=TRUE),
  total_transcript_observations=sum(transcripts),
  total_FL_support=sum(total_FL_support,na.rm=TRUE)),by=group]
save_table(group_cat,'group_structural_categories.tsv')
save_table(cohort_cat,'cohort_structural_categories.tsv')
save_table(group_yield,'group_yield_summary.tsv')

# --- 5. Figures --------------------------------------------------------------
# 01: 274 individual samples, ordered by group and within group by yield.
setorder(samples,group,transcripts)
ordered_ids <- samples$sample_id
prop <- copy(cats)
prop[,proportion:=transcripts/sum(transcripts),by=sample_id]
prop[,sample_id:=factor(sample_id,levels=ordered_ids)]
p1 <- ggplot(prop,aes(sample_id,proportion,fill=category))+
  geom_col(width=1)+scale_category()+scale_y_continuous(labels=percent_format())+
  labs(title='SQANTI3 structural categories across 274 samples',
       subtitle='Rescued transcript models before IsoMatch; samples ordered within disease group',
       x='Individual samples',y='Transcript models (%)')+
  facet_grid(~group,scales='free_x',space='free_x')+plot_theme+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.spacing.x=unit(0.2,'lines'))
save_plot(p1,'Figure_01_sample_category_proportions',width=17,height=7)
# 02: Requested group-wise stacked bars (matches example palette and style).
p2 <- ggplot(group_cat,aes(group,transcripts,fill=category))+
  geom_col(position='fill',width=.78)+scale_category()+
  scale_y_continuous(labels=percent_format(accuracy=1),expand=expansion(mult=c(0,.02)))+
  labs(title='Structural category composition by disease group',
       subtitle='Proportion of rescued sample-level transcript observations (not unique cohort isoforms)',
       x=NULL,y='Transcript models (%)')+plot_theme
save_plot(p2,'Figure_02_group_category_100pct',width=11,height=7)
# 03: Group absolute counts (unequal sample sizes are NOT normalized).
p3 <- ggplot(group_cat,aes(group,transcripts,fill=category))+
  geom_col(width=.78)+scale_category()+scale_y_continuous(labels=label_comma())+
  labs(title='Number of transcript observations by disease group',
       subtitle='Absolute totals; group sizes differ, so compare with proportional composition and per-sample yield',
       x=NULL,y='Rescued transcript observations')+plot_theme
save_plot(p3,'Figure_03_group_category_absolute',width=11,height=7)
# 04: Whole-cohort structural composition, count vs FL support.
cohort_long <- melt(copy(cohort_cat)[, transcripts := as.numeric(transcripts)],id.vars='category',
  measure.vars=c('transcripts','FL_support'),variable.name='measure',value.name='count')
cohort_long <- cohort_long[is.finite(count)]
cohort_long[,measure:=factor(measure,levels=c('transcripts','FL_support'),
                             labels=c('Transcript observations','Aggregated FL-read support'))]
p4 <- ggplot(cohort_long,aes(category,count,fill=category))+
  geom_col()+facet_wrap(~measure,scales='free_y',ncol=1)+scale_category()+
  scale_y_continuous(labels=label_comma())+
  labs(title='Cohort-wide SQANTI3 structural categories',
       subtitle='Summed over 274 samples; repeated transcript structures across samples are counted repeatedly',
       x=NULL,y='Count')+plot_theme+theme(axis.text.x=element_text(angle=35,hjust=1))
save_plot(p4,'Figure_04_cohort_category_counts',width=13,height=10)
# 04A / 04B: Two pie charts from the SAME cohort summary as the cohort-wide summary.
# The cohort_cat table is created in section 4 and has columns:
# category, transcripts, FL_support. Do not use cohort_category.
# Counts sum sample-level SQANTI3 observations (not unique isoforms).
make_category_pie <- function(dt, value_col, title, subtitle) {
  d <- data.table(
    category = as.character(dt$category),
    count = as.numeric(dt[[value_col]])
  )
  d <- d[is.finite(count) & count > 0]
  if (!nrow(d) || sum(d$count) <= 0) {
    warning('Skipping pie chart: no positive values for ', value_col)
    return(NULL)
  }
  d[, percent := 100 * count / sum(count)]
  d[, category := factor(category, levels=CATEGORY_LEVELS)]
  setorder(d, category)
  # Use identical category vectors for legend breaks and labels.
  # Only include categories actually observed with positive counts.
  present_categories <- CATEGORY_LEVELS[CATEGORY_LEVELS %in% as.character(d$category)]
  present_pct <- d$percent[match(present_categories, as.character(d$category))]
  legend_labels <- paste0(present_categories, ' (', sprintf('%.1f%%', present_pct), ')')
  stopifnot(length(present_categories) == length(legend_labels))
  ggplot(d, aes(x=1, y=count, fill=category)) +
    geom_col(width=1, color='white', linewidth=0.25) +
    coord_polar(theta='y') +
    scale_fill_manual(values=CATEGORY_COLORS,
                      breaks=present_categories,
                      labels=legend_labels,
                      drop=TRUE, name='Structural category') +
    xlim(0.5, 1.5) +
    labs(title=title, subtitle=subtitle) +
    theme_void(base_size=12) +
    theme(plot.title=element_text(face='bold',size=14),
          plot.subtitle=element_text(size=10),
          legend.position='right')
}
p4a <- make_category_pie(
  cohort_cat, 'transcripts',
  'Cohort SQANTI3 transcript observations',
  '274 samples; repeated transcript structures counted once per sample observation'
)
if (!is.null(p4a)) save_plot(p4a, 'Figure_04A_cohort_observed_pie', width=11, height=7)
p4b <- make_category_pie(
  cohort_cat, 'FL_support',
  'Cohort SQANTI3 FL-read support',
  '274 samples; structural categories weighted by aggregated FL support'
)
if (!is.null(p4b)) save_plot(p4b, 'Figure_04B_cohort_FL_pie', width=11, height=7)
pie_summary <- melt(
  copy(cohort_cat)[, transcripts := as.numeric(transcripts)],
  id.vars='category', measure.vars=c('transcripts','FL_support'),
  variable.name='metric', value.name='count'
)
pie_summary[, percent := if (sum(count,na.rm=TRUE)>0)
  100 * count / sum(count,na.rm=TRUE) else NA_real_, by=metric]
save_table(pie_summary, 'Figure_04_pie_category_percentages.tsv')
# 05: Per-sample yield grouped by disease; avoids confounding group size.
sample_long <- melt(samples,id.vars=c('sample_id','group','sex'),
  measure.vars=c('transcripts','total_FL_support'),variable.name='measure',value.name='value')
sample_long[,measure:=factor(measure,levels=c('transcripts','total_FL_support'),
  labels=c('Transcript observations','Aggregated FL-read support'))]
p5 <- ggplot(sample_long,aes(group,value))+
  geom_boxplot(width=.65,fill='#AFC8DF',outlier.shape=NA)+
  geom_jitter(width=.15,height=0,alpha=.45,size=.9,color='#29496C')+
  facet_wrap(~measure,scales='free_y',ncol=1)+scale_y_continuous(labels=label_comma())+
  labs(title='Sample-level transcript yield by disease group',
       subtitle='Each point is one individual; boxplots show group medians and quartiles',
       x=NULL,y='Count per sample')+plot_theme
save_plot(p5,'Figure_05_group_sample_yield',width=12,height=10)
# 06: FL support weighted group proportions, only if FL is available.
if(any(is.finite(group_cat$FL_support)) && sum(group_cat$FL_support,na.rm=TRUE)>0) {
  p6 <- ggplot(group_cat[is.finite(FL_support)],aes(group,FL_support,fill=category))+
    geom_col(position='fill',width=.78)+scale_category()+
    scale_y_continuous(labels=percent_format(accuracy=1),expand=expansion(mult=c(0,.02)))+
    labs(title='FL-read support by disease group and structural class',
         subtitle='Each class weighted by aggregated FL support rather than number of transcript models',
         x=NULL,y='FL-read support (%)')+plot_theme
  save_plot(p6,'Figure_06_group_FL_weighted_category_100pct',width=11,height=7)
}
# 07: Sample-level category distributions by disease group.
prop[,group:=factor(group,levels=GROUP_LABELS)]
p7 <- ggplot(prop,aes(group,proportion,fill=category))+
  geom_boxplot(outlier.alpha=.18,width=.7)+
  facet_wrap(~category,ncol=3,scales='free_y')+
  scale_category()+scale_y_continuous(labels=percent_format())+
  labs(title='Between-sample variation in structural categories',
       subtitle='Within-sample category percentages; boxplots summarize individuals in each group',
       x=NULL,y='Within-sample fraction')+plot_theme+
  theme(legend.position='none',axis.text.x=element_text(angle=35,hjust=1))
save_plot(p7,'Figure_07_group_sample_category_distributions',width=14,height=12)
# 08: Length summaries by group/category, when SQANTI3 length exists.
if(nrow(lengths)) {
  p8 <- ggplot(lengths,aes(group,median_length,fill=category))+
    geom_boxplot(outlier.alpha=.15)+facet_wrap(~category,ncol=3,scales='free_y')+
    scale_category()+scale_y_continuous(labels=label_comma())+
    labs(title='Transcript length by disease group and SQANTI3 class',
         subtitle='Distribution of per-sample median lengths, not individual transcript lengths',
         x=NULL,y='Median transcript length (nt)')+plot_theme+
    theme(legend.position='none',axis.text.x=element_text(angle=35,hjust=1))
  save_plot(p8,'Figure_08_group_transcript_lengths',width=14,height=12)
}
# 09: Quality indicators by group, only for fields actually present.
if(nrow(quality)) {
  qplot_data <- quality[is.finite(fraction)&fraction>=0&fraction<=1]
  if(nrow(qplot_data)) {
    p9 <- ggplot(qplot_data,aes(group,fraction))+
      geom_boxplot(fill='#AFC8DF',outlier.alpha=.2)+
      facet_wrap(~metric,scales='free_y')+scale_y_continuous(labels=percent_format())+
      labs(title='Available SQANTI3 quality indicators by disease group',
           subtitle='Only metrics present in the rescued classification files are plotted',
           x=NULL,y='Fraction of transcript models')+plot_theme+
      theme(axis.text.x=element_text(angle=35,hjust=1))
    save_plot(p9,'Figure_09_group_quality_indicators',width=13,height=8)
  }
}
# 10: Optional sex-by-group structural composition.
if(MAKE_SEX_STRATIFIED) {
  sex_cat <- cats[,.(transcripts=sum(transcripts)),by=.(group,sex,category)]
  save_table(sex_cat,'group_sex_structural_categories.tsv')
  p10 <- ggplot(sex_cat,aes(group,transcripts,fill=category))+
    geom_col(position='fill',width=.78)+facet_wrap(~sex,ncol=1)+scale_category()+
    scale_y_continuous(labels=percent_format())+
    labs(title='SQANTI3 structural categories by group and sex',
         subtitle='Pooled transcript observations within each disease group and sex',
         x=NULL,y='Transcript models (%)')+plot_theme
  save_plot(p10,'Figure_10_group_sex_category_100pct',width=11,height=10)
}

# --- 6. Inventory of upstream SQANTI3 filtering stages ----------------------
# Do not infer prefilter/filtered/rescued counts from rescued-only files.
all_class <- list.files(PREMERGE_DIR,pattern='classification\\.txt(\\.gz)?$',
                        recursive=TRUE,full.names=TRUE)
all_class <- all_class[!grepl('sqanti3_with_sr|with_sr',all_class,ignore.case=TRUE)]
inv <- data.table(sample_id=manifest_ids)
for(stage in c('initial','rules_filtered','rescued')) {
  suffix <- switch(stage,initial='_classification.txt',
                   rules_filtered='_RulesFilter_classification.txt',
                   rescued='_rescued_classification.txt')
  matches <- all_class[grepl(paste0(gsub('\\.','\\\\.',suffix),'(\\.gz)?$'),all_class)]
  if(stage=='initial') matches <- matches[!grepl('RulesFilter|rescued',basename(matches),ignore.case=TRUE)]
  found <- sub(paste0(gsub('\\.','\\\\.',suffix),'(\\.gz)?$'),'',basename(matches))
  inv[,(paste0('has_',stage)):=sample_id %in% found]
}
save_table(inv,'SQANTI3_filter_stage_file_inventory.tsv')
message('DONE. Figures and tables saved to: ',OUT)
message('All results above are PRE-MERGE rescued sample-level SQANTI3 observations.')
message('Stage-specific filtering counts require actual initial and filtered classification files.')
