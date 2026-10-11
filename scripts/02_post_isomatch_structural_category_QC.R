# ============================================================================
# CWOW | POST-ISOMATCH MODERATE-CONFIDENCE STRUCTURAL CATEGORY QC | 274 samples
# Run interactively in Posit Workbench. No command-line arguments.
#
# PURPOSE
#   Describe the EXACT transcript IDs retained in the moderate-confidence
#   count matrix used for DTE/DTU. Recover source SQANTI3 structural classes
#   using IsoMatch tracking + rescued sample-level SQANTI3 classifications.
#   No post-merge SQANTI3 reclassification is implied or performed.
#
# CLASSIFICATION RULES
#   All matched source classifications agree -> that SQANTI3 category.
#   Different matched source categories -> Mixed (reported separately).
#   No matched source classification -> Unresolved (reported separately).
#   Partial matches are flagged; category assignments with unmatched source
#   mappings should be interpreted cautiously. No majority-vote guessing.
#   Source transcript observations and unique merged models are NOT equivalent.
#
# INPUTS
#   1. Moderate-confidence DTE/DTU transcript count matrix (.tsv.gz)
#   2. IsoMatch merged transcript tracking file (.tsv.gz)
#   3. Original 274 rescued SQANTI3 classifications
#   4. Cohort metadata (sample ID, disease group, sex)
#
# OUTPUT DIRECTORY: OUT
# TABLES (.tsv)
#   retained_transcript_structural_categories.tsv   One row per retained ISOMT
#   source_mapping_diagnostics.tsv                 Source/ID match diagnostics
#   cohort_category_composition.tsv                Unique retained models by class
#   sample_category_composition.tsv                Per-sample detected models + FL
#   group_category_composition.tsv                 Disease-group detected models + FL
#   group_sex_category_composition.tsv             Disease-by-sex detected models + FL
#   sample_metadata_used.tsv                       Count matrix metadata alignment
#   category_assignment_overview.tsv               Category assignment QC counts
# FIGURES (.png and .pdf)
#   Figure_01_cohort_unique_category_pie           Unique retained ISOMT composition
#   Figure_02_cohort_unique_category_counts        Unique retained ISOMT counts
#   Figure_03_group_category_FL_100pct             FL-weighted disease composition
#   Figure_04_group_category_detected_100pct       Detected transcript observations
#   Figure_05_sample_category_FL_100pct            Sample FL-weighted composition
#   Figure_06_group_sex_category_FL_100pct         Disease x sex FL composition
#   Figure_07_sample_category_distributions        Per-sample class fractions
#   Figure_08_category_assignment_QC               Agreement/mapping outcomes
#
# IMPORTANT
#   These are TRACKING-DERIVED source SQANTI3 labels, NOT a de novo SQANTI3
#   classification of merged transcript structures. Inspect diagnostic tables.
#   Group-level detected transcript counts count a transcript once PER SAMPLE
#   where FL > 0, while cohort unique counts count each ISOMT exactly once.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# --- User settings ----------------------------------------------------------
ROOT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/'
BASE <- file.path(ROOT, 'results/isomatch_outputs')
COUNTS <- file.path(BASE, '04_confidence_filter_rescuedLR_274/moderate/CWOW_274_rescuedLR.confidence_filtered.counts.tsv.gz')
TRACK <- file.path(BASE, '02_out_CWOWSamplesCohort/CWOWSamplesCohort.track.tsv.gz')
PREMERGE_DIR <- file.path(ROOT, 'scripts/grandcanyon/results')
METADATA <- file.path(ROOT, 'scripts/CWOW_full_long_read_cohort_n276.csv')
OUT <- file.path(ROOT, 'results/postmerge_moderate_structural_QC')
EXPECTED_SAMPLES <- 274L
CHUNK <- 150000L
DPI <- 300
CATEGORY_LEVELS <- c('FSM','ISM','NIC','NNC','Fusion','Genic','Genic intron',
                     'Antisense','Intergenic','More junctions','Other','Mixed','Unresolved')
CATEGORY_COLORS <- c('FSM'='#4F83B6','ISM'='#E99545','NIC'='#F5C342','NNC'='#9D9FA1',
                     'Fusion'='#75B1D0','Genic'='#79AA5B','Genic intron'='#B99A73',
                     'Antisense'='#29496C','Intergenic'='#963B20',
                     'More junctions'='#C65C91','Other'='#B8B8B8',
                     'Mixed'='#7651A8','Unresolved'='#DADADA')
GROUP_LEVELS <- c('CONTROL','AD','LBD_S','LBD_AS','LBD_ATS')
GROUP_LABELS <- c('Control','AD(AT)','LBD(S)','LBD(AS)','LBD(ATS)')
for (p in c(COUNTS, TRACK, METADATA)) if (!file.exists(p)) stop('Missing input: ',p)
if (!dir.exists(PREMERGE_DIR)) stop('Missing classification directory: ',PREMERGE_DIR)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
save_tsv <- function(d, file) fwrite(d,file.path(OUT,file),sep='\t',na='NA')
save_plot <- function(p, stem, w=11, h=7) {
  ggsave(file.path(OUT,paste0(stem,'.png')),p,width=w,height=h,dpi=DPI,bg='white',limitsize=FALSE)
  ggsave(file.path(OUT,paste0(stem,'.pdf')),p,width=w,height=h,bg='white',limitsize=FALSE)
  message('Saved: ',stem)
}
normalize_class <- function(x) {
  x <- tolower(trimws(as.character(x)))
  fcase(x %in% c('full-splice_match','full-splice match','fsm'),'FSM',
        x %in% c('incomplete-splice_match','incomplete-splice match','ism'),'ISM',
        x %in% c('novel_in_catalog','novel in catalog','nic'),'NIC',
        x %in% c('novel_not_in_catalog','novel not in catalog','nnc'),'NNC',
        grepl('fusion',x),'Fusion',
        x %in% c('genic','genic_genomic','genic-genomic'),'Genic',
        x %in% c('genic_intron','genic-intron'),'Genic intron',
        grepl('antisense',x),'Antisense',grepl('intergenic',x),'Intergenic',
        grepl('more.*junction',x),'More junctions',default='Other')
}
plot_theme <- theme_classic(base_size=12)+theme(plot.title=element_text(face='bold'),
  axis.text.x=element_text(color='black'),legend.position='right')
category_scale <- function() scale_fill_manual(values=CATEGORY_COLORS,drop=FALSE,name='Structural category')

# --- 1. Read exact retained transcript catalog and metadata -----------------
message('Reading retained moderate-confidence count matrix...')
counts <- fread(COUNTS,showProgress=TRUE)
setnames(counts,1L,'transcript_id')
if (anyNA(counts$transcript_id) || anyDuplicated(counts$transcript_id))
  stop('Retained count matrix has missing or duplicated transcript IDs.')
sample_ids <- setdiff(names(counts),'transcript_id')
if (length(sample_ids)!=EXPECTED_SAMPLES)
  warning('Expected ',EXPECTED_SAMPLES,' samples; count matrix has ',length(sample_ids))
meta <- fread(METADATA)
if (!all(c('Study-Subject-ID','group','Gender') %in% names(meta)))
  stop('Metadata must contain Study-Subject-ID, group, Gender')
meta[,sample_id:=trimws(as.character(`Study-Subject-ID`))]
meta[!grepl('^CWOW_',sample_id),sample_id:=paste0('CWOW_',sample_id)]
meta[,group:=toupper(gsub('[- ]','_',trimws(as.character(group))))]
meta[,sex:=toupper(trimws(as.character(Gender)))]
meta[sex %in% c('F','FEMALE','XX'),sex:='XX']
meta[sex %in% c('M','MALE','XY'),sex:='XY']
meta <- meta[sample_id %in% sample_ids,.(sample_id,group,sex)]
if (anyDuplicated(meta$sample_id) || !setequal(meta$sample_id,sample_ids))
  stop('Sample IDs in metadata do not uniquely match count matrix column names.')
if (any(!meta$group %in% GROUP_LEVELS) || any(!meta$sex %in% c('XX','XY')))
  stop('Unexpected disease group or sex in metadata.')
meta[,group:=factor(group,levels=GROUP_LEVELS,labels=GROUP_LABELS)]
meta[,sex:=factor(sex,levels=c('XX','XY'))]
setkey(meta,sample_id)
save_tsv(meta,'sample_metadata_used.tsv')
retained_ids <- counts$transcript_id
message('Retained transcripts: ',format(length(retained_ids),big.mark=','))

# --- 2. Identify rescued sample SQANTI3 files -------------------------------
files <- list.files(PREMERGE_DIR,pattern='_rescued_classification\\.txt(\\.gz)?$',
                    recursive=TRUE,full.names=TRUE)
files <- files[!grepl('sqanti3_with_sr|with_sr',files,ignore.case=TRUE)]
file_ids <- sub('_rescued_classification\\.txt(\\.gz)?$','',basename(files))
file_index <- data.table(sample_id=file_ids,path=files)[sample_id %in% sample_ids]
if (anyDuplicated(file_index$sample_id)) stop('Duplicate rescued SQANTI3 files for sample(s).')
if (!setequal(file_index$sample_id,sample_ids))
  stop('Missing rescued SQANTI3 files: ',paste(setdiff(sample_ids,file_index$sample_id),collapse=', '))
setkey(file_index,sample_id)

# --- 3. Stream IsoMatch tracking; keep retained ISOMT source mappings -------
message('Streaming IsoMatch tracking file...')
con <- gzfile(TRACK,'rt')
header <- strsplit(readLines(con,n=1L),'\t',fixed=TRUE)[[1L]]
needed <- c('merged_tx_id','src_tx_id','src_file_name')
if (!all(needed %in% header)) stop('Tracking columns missing: ',paste(setdiff(needed,header),collapse=', '))
pos <- match(needed,header)
parts <- list(); k <- 0L; n_processed <- 0L
repeat {
  lines <- readLines(con,n=CHUNK,warn=FALSE)
  if (!length(lines)) break
  fields <- strsplit(lines,'\t',fixed=TRUE)
  valid <- lengths(fields)>=max(pos)
  if (any(!valid)) warning('Malformed tracking rows skipped in this chunk: ',sum(!valid))
  fields <- fields[valid]
  merged <- vapply(fields,`[`,character(1),pos[1])
  hits <- which(merged %chin% retained_ids)
  if (length(hits)) {
    parts[[k+1L]] <- data.table(transcript_id=merged[hits],
      src_tx_id=vapply(fields[hits],`[`,character(1),pos[2]),
      src_file_name=vapply(fields[hits],`[`,character(1),pos[3]))
    k <- k+1L
  }
  n_processed <- n_processed+length(lines)
  if (n_processed %% 3000000L < length(lines)) message('Tracking lines processed: ',format(n_processed,big.mark=','))
}
close(con)
if (!length(parts)) stop('No retained ISOMT IDs found in tracking file. Check inputs.')
track <- unique(rbindlist(parts,use.names=TRUE))
rm(parts); gc(verbose=FALSE)

# Match source filenames to sample IDs using full known sample names, not a
# guessed regular expression. Sample IDs are matched as basename prefixes.
source_files <- unique(track$src_file_name)
source_index <- data.table(src_file_name=source_files,sample_id=NA_character_)
known_ids <- sample_ids[order(nchar(sample_ids),decreasing=TRUE)]
for (id in known_ids) {
  base <- basename(source_index$src_file_name)
  hit <- is.na(source_index$sample_id) & (base==id | startsWith(base,paste0(id,'.')) |
           startsWith(base,paste0(id,'_')) | startsWith(base,paste0(id,'-')))
  source_index[hit,sample_id:=id]
}
track <- merge(track,source_index,by='src_file_name',all.x=TRUE,sort=FALSE)
message('Tracking source mappings: ',format(nrow(track),big.mark=','),
        '; unresolved source filenames: ',sum(is.na(track$sample_id)))

# --- 4. Match source transcript IDs to rescued SQANTI3 classifications ------
# Process one sample at a time; retain only transcript IDs present in tracking.
track[,category:=NA_character_]
for (i in seq_len(nrow(file_index))) {
  id <- file_index$sample_id[i]
  idx <- which(!is.na(track$sample_id) & track$sample_id==id)
  if (!length(idx)) next
  path <- file_index$path[i]
  cols <- names(fread(path,nrows=0,showProgress=FALSE))
  if (!all(c('isoform','structural_category') %in% cols)) stop('Missing SQANTI3 columns: ',path)
  d <- fread(path,select=c('isoform','structural_category'),showProgress=FALSE)
  if (anyDuplicated(d$isoform)) stop('Duplicated SQANTI3 isoform IDs in: ',path)
  wanted <- unique(track$src_tx_id[idx])
  d <- d[isoform %chin% wanted]
  d[,category:=normalize_class(structural_category)]
  track[idx,category:=d$category[match(track$src_tx_id[idx],d$isoform)]]
  if (i %% 20L==0L || i==nrow(file_index)) message('Matched SQANTI3 files: ',i,'/',nrow(file_index))
}

# --- 5. One category per retained ISOMT, with transparent ambiguity ----------
# Aggregate exact source mappings. A single merged ISOMT may have multiple
# distinct SQANTI3 categories; never force a majority category.
assignment <- track[,{
  observed <- unique(category[!is.na(category)])
  n_matched <- sum(!is.na(category))
  .(n_source_mappings=.N,n_matched_source_mappings=n_matched,
    n_unmatched_source_mappings=.N-n_matched,
    n_distinct_source_categories=length(observed),
    source_categories=paste(sort(observed),collapse=';'),
    category=if (!length(observed)) 'Unresolved' else if (length(observed)>1L) 'Mixed' else observed,
    assignment_status=if (!length(observed)) 'Unresolved' else if (length(observed)>1L) 'Mixed' else if (n_matched<.N) 'Partial match' else 'Concordant')
},by=transcript_id]
assignment <- merge(data.table(transcript_id=retained_ids),assignment,by='transcript_id',all.x=TRUE,sort=FALSE)
assignment[is.na(category),`:=`(category='Unresolved',assignment_status='No tracking record',
                                n_source_mappings=0L,n_matched_source_mappings=0L,
                                n_unmatched_source_mappings=0L,n_distinct_source_categories=0L,
                                source_categories='')]
assignment[,category:=factor(category,levels=CATEGORY_LEVELS)]
save_tsv(assignment,'retained_transcript_structural_categories.tsv')
save_tsv(assignment[,.(n_transcripts=.N),by=.(assignment_status,category)],'category_assignment_overview.tsv')
save_tsv(data.table(metric=c('retained_ISOMT','tracking_source_mappings',
  'unmatched_source_file_mappings','unmatched_SQANTI3_source_mappings',
  'retained_ISOMT_without_tracking','retained_ISOMT_mixed','retained_ISOMT_partial'),
  value=c(nrow(assignment),nrow(track),sum(is.na(track$sample_id)),
          sum(is.na(track$category)),sum(assignment$assignment_status=='No tracking record'),
          sum(assignment$assignment_status=='Mixed'),sum(assignment$assignment_status=='Partial match'))),
  'source_mapping_diagnostics.tsv')
if (mean(!is.na(track$category))<0.5)
  warning('Fewer than 50% of source tracking mappings matched SQANTI3. Inspect source_mapping_diagnostics.tsv; check source IDs and file naming.')

# --- 6. Unique cohort composition and count-weighted sample summaries -------
setkey(assignment,transcript_id)
setkey(counts,transcript_id)
counts <- assignment[counts,on='transcript_id']
if (anyNA(counts$category)) stop('Unassigned retained transcript category encountered.')
counts[,category:=factor(as.character(category),levels=CATEGORY_LEVELS)]
cohort <- counts[,.(n_unique_transcripts=.N),by=category]
cohort[,pct_unique:=100*n_unique_transcripts/sum(n_unique_transcripts)]
save_tsv(cohort,'cohort_category_composition.tsv')

# One row per sample and category. FL counts are not transcript counts.
# A detected transcript observation requires FL > 0 in that sample.
sample_parts <- vector('list',length(sample_ids))
for (i in seq_along(sample_ids)) {
  id <- sample_ids[i]
  x <- suppressWarnings(as.numeric(counts[[id]]))
  if (any(!is.finite(x)) || any(x<0)) stop('Nonfinite or negative FL count: ',id)
  sample_parts[[i]] <- data.table(category=counts$category,FL=x)[,
    .(detected_transcripts=sum(FL>0),FL_support=sum(FL)),by=category][,sample_id:=id][]
  if (i %% 25L==0L || i==length(sample_ids)) message('Aggregated count columns: ',i,'/',length(sample_ids))
}
sample_cat <- rbindlist(sample_parts)
sample_cat <- merge(sample_cat,meta,by='sample_id',all.x=TRUE,sort=FALSE)
sample_cat[,pct_detected:=if(sum(detected_transcripts)>0) 100*detected_transcripts/sum(detected_transcripts) else NA_real_,by=sample_id]
sample_cat[,pct_FL:=if(sum(FL_support)>0) 100*FL_support/sum(FL_support) else NA_real_,by=sample_id]
save_tsv(sample_cat,'sample_category_composition.tsv')
group_cat <- sample_cat[,.(detected_transcripts=sum(detected_transcripts),FL_support=sum(FL_support)),by=.(group,category)]
group_cat[,pct_detected:=100*detected_transcripts/sum(detected_transcripts),by=group]
group_cat[,pct_FL:=100*FL_support/sum(FL_support),by=group]
save_tsv(group_cat,'group_category_composition.tsv')
sex_cat <- sample_cat[,.(detected_transcripts=sum(detected_transcripts),FL_support=sum(FL_support)),by=.(group,sex,category)]
sex_cat[,pct_detected:=100*detected_transcripts/sum(detected_transcripts),by=.(group,sex)]
sex_cat[,pct_FL:=100*FL_support/sum(FL_support),by=.(group,sex)]
save_tsv(sex_cat,'group_sex_category_composition.tsv')

# --- 7. Figures (same SQANTI3 palette as pre-merge QC) ----------------------
make_pie <- function(d) {
  d <- copy(d[n_unique_transcripts>0])
  d[,category:=factor(as.character(category),levels=CATEGORY_LEVELS)]
  present <- as.character(d$category)
  labels <- paste0(present,' (',sprintf('%.1f%%',d$pct_unique),')')
  ggplot(d,aes(x=1,y=n_unique_transcripts,fill=category))+
    geom_col(width=1,color='white',linewidth=.25)+coord_polar(theta='y')+
    scale_fill_manual(values=CATEGORY_COLORS,breaks=present,labels=labels,drop=TRUE)+
    xlim(.5,1.5)+labs(title='Moderate-confidence merged transcript catalog',
      subtitle='Unique retained ISOMT models; tracking-derived SQANTI3 categories')+
    theme_void(base_size=12)+theme(plot.title=element_text(face='bold'))
}
save_plot(make_pie(cohort),'Figure_01_cohort_unique_category_pie')
p2 <- ggplot(cohort,aes(category,n_unique_transcripts,fill=category))+
  geom_col()+category_scale()+scale_y_continuous(labels=label_comma())+
  labs(title='Structural composition of retained merged transcript models',
       subtitle='Each moderate-confidence ISOMT counted once',x=NULL,y='Unique transcripts')+
  plot_theme+theme(axis.text.x=element_text(angle=35,hjust=1))
save_plot(p2,'Figure_02_cohort_unique_category_counts',12,7)
p3 <- ggplot(group_cat,aes(group,FL_support,fill=category))+
  geom_col(position='fill',width=.78)+category_scale()+
  scale_y_continuous(labels=percent_format())+
  labs(title='Post-merge structural categories by disease group',
       subtitle='Aggregated FL counts from the moderate-confidence DTE/DTU input matrix',
       x=NULL,y='FL-read support (%)')+plot_theme
save_plot(p3,'Figure_03_group_category_FL_100pct')
p4 <- ggplot(group_cat,aes(group,detected_transcripts,fill=category))+
  geom_col(position='fill',width=.78)+category_scale()+
  scale_y_continuous(labels=percent_format())+
  labs(title='Detected merged transcripts by disease group',
       subtitle='Each ISOMT counted once per sample where FL > 0; not unique within group',
       x=NULL,y='Detected transcript observations (%)')+plot_theme
save_plot(p4,'Figure_04_group_category_detected_100pct')
sample_cat[,sample_order:=factor(sample_id,levels=meta[order(group,sample_id),sample_id])]
p5 <- ggplot(sample_cat,aes(sample_order,FL_support,fill=category))+
  geom_col(position='fill',width=1)+category_scale()+facet_grid(~group,scales='free_x',space='free_x')+
  scale_y_continuous(labels=percent_format())+
  labs(title='Moderate-confidence structural composition across samples',
       subtitle='FL-weighted; each bar is one individual',x='Individual samples',y='FL-read support (%)')+
  plot_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
save_plot(p5,'Figure_05_sample_category_FL_100pct',17,7)
p6 <- ggplot(sex_cat,aes(group,FL_support,fill=category))+
  geom_col(position='fill',width=.78)+facet_wrap(~sex,ncol=1)+category_scale()+
  scale_y_continuous(labels=percent_format())+
  labs(title='Post-merge structural categories by disease group and sex',
       subtitle='Moderate-confidence transcripts; aggregated FL-read support',
       x=NULL,y='FL-read support (%)')+plot_theme
save_plot(p6,'Figure_06_group_sex_category_FL_100pct',11,10)
p7 <- ggplot(sample_cat,aes(group,pct_FL/100,fill=category))+
  geom_boxplot(outlier.alpha=.15)+facet_wrap(~category,ncol=3,scales='free_y')+
  category_scale()+scale_y_continuous(labels=percent_format())+
  labs(title='Between-sample variation in structural category support',
       subtitle='Within-sample FL fractions; one observation per individual',
       x=NULL,y='Within-sample FL fraction')+plot_theme+
  theme(legend.position='none',axis.text.x=element_text(angle=35,hjust=1))
save_plot(p7,'Figure_07_sample_category_distributions',14,12)
qc <- assignment[,.(n_transcripts=.N),by=assignment_status]
p8 <- ggplot(qc,aes(reorder(assignment_status,n_transcripts),n_transcripts,fill=assignment_status))+
  geom_col(show.legend=FALSE)+coord_flip()+scale_y_continuous(labels=label_comma())+
  labs(title='Source classification assignment quality',
       subtitle='Concordance and completeness of IsoMatch-to-SQANTI3 mapping',
       x=NULL,y='Retained merged transcripts')+theme_classic(base_size=12)
save_plot(p8,'Figure_08_category_assignment_QC')
message('DONE: ',OUT)
message('Review category_assignment_overview.tsv and source_mapping_diagnostics.tsv before interpreting category plots.')
