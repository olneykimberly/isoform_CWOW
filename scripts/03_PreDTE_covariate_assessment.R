# ============================================================================
# CWOW | Pre-DTE technical QC, RIN, PCA, CCA and variance partition
# Source directly in Posit Workbench. Read-only; does not alter DTE inputs.
# INPUTS: moderate-confidence FL count matrix; CWOW cohort metadata.
# OUTPUT: results/QC_technical_RIN_normalization_CCA
# FIGURES: RIN/group violin; RIN associations; normalization diagnostics with
# outlier IDs; PCA by group, sex, SMRT cell and quantitative traits; group-wise
# Kruskal-Wallis violins; CCA of expression PC1-10 versus clinical variables;
# per-transcript variance partition (optional variancePartition package).
# SMRT cell RGIDs are parsed from BOTH dataFiles in each sample JSON.
# Missing JSONs are explicitly reported; sample IDs are validated before variance partition.
# Shared SMRT/run IDs are exploratory technical covariates; never assign one
# of two SMRT cells as a unique sample-level batch without validation.
# Pathology variables are exploratory correlates, not automatic DTE covariates.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(DESeq2); library(matrixStats); library(jsonlite); library(readxl)
})
BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
COUNTS <- file.path(BASE,'04_confidence_filter_rescuedLR_274','moderate','CWOW_274_rescuedLR.confidence_filtered.counts.tsv.gz')
METADATA <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/CWOW_full_long_read_cohort_n276.csv'
OUT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/QC_technical_RIN_normalization_CCA/'
RNA_SAMPLES_XLSX <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/metadata/CWOW_Original_RNA_Samples.xlsx'
JSON_DIR <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/output/jsons'
MIN_FL <- 5L; MIN_SAMPLES <- 5L; TOP_VARIABLE <- 2000L
VARPART_TRANSCRIPTS <- 10000L; RUN_VARPART <- TRUE
GROUP_LEVELS <- c('CONTROL','AD','LBD_S','LBD_AS','LBD_ATS')
GROUP_LABELS <- c('Control','AD(AT)','LBD(S)','LBD(AS)','LBD(ATS)')
GROUP_COLORS <- c('Control'='#377EB8','AD(AT)'='#E41A1C','LBD(S)'='#6A3D9A','LBD(AS)'='#9E9AC8','LBD(ATS)'='#636363')
SEX_COLORS <- c('Female'='#8A5DB7','Male'='#E99545')
for(p in c(COUNTS,METADATA,RNA_SAMPLES_XLSX)) if(!file.exists(p)) stop('Missing input: ',p)
if(!dir.exists(JSON_DIR)) stop('Missing JSON directory: ',JSON_DIR)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write_tab <- function(x,n) fwrite(x,file.path(OUT,n),sep='\t',na='NA')
save_plot <- function(p,n,w=9,h=6) {
  ggsave(file.path(OUT,paste0(n,'.png')),p,width=w,height=h,dpi=300,bg='white',limitsize=FALSE)
  ggsave(file.path(OUT,paste0(n,'.pdf')),p,width=w,height=h,bg='white',limitsize=FALSE)
  message('Saved ',n)
}
base_theme <- theme_bw(base_size=12)+theme(plot.title=element_text(face='bold'),axis.text.x=element_text(angle=25,hjust=1))
scale_group_fill <- function() scale_fill_manual(values=GROUP_COLORS,drop=FALSE)
scale_group_color <- function() scale_color_manual(values=GROUP_COLORS,drop=FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
# Metadata: retain original columns, normalize labels and construct plot ID.
md <- fread(METADATA,check.names=FALSE)
# Preserve every original metadata field unchanged in the augmented cohort export.
md_original <- copy(md)
if(!'Study-Subject-ID'%in%names(md)) stop('Metadata lacks Study-Subject-ID')
md[,sample_id:=paste0('CWOW_',sub('^CWOW_','',trimws(as.character(`Study-Subject-ID`))))]
if(anyDuplicated(md$sample_id)) stop('Duplicate metadata sample IDs')
required <- c('group','Gender','Tapestation_DIN_or_RIN','Brain.wt','PMI','Braak.NFT','Thal.amyloid','Cing.LB','Age','SMRTcell_ID')
if(length(setdiff(required,names(md)))) stop('Missing metadata fields: ',paste(setdiff(required,names(md)),collapse=', '))
md[,group_code:=toupper(gsub('-','_',trimws(as.character(group))))]
if(any(!md$group_code%in%GROUP_LEVELS,na.rm=TRUE)) stop('Unexpected disease-group labels')
md[,group:=factor(group_code,levels=GROUP_LEVELS,labels=GROUP_LABELS)]
md[,sex:=fcase(toupper(trimws(as.character(Gender)))%in%c('F','FEMALE','XX'),'Female',toupper(trimws(as.character(Gender)))%in%c('M','MALE','XY'),'Male',default=NA_character_)]
md[,sex:=factor(sex,levels=c('Female','Male'))]
for(pair in list(c('Brain.wt','Brain_weight'),c('Braak.NFT','Braak'),c('Thal.amyloid','Thal'),c('Cing.LB','CingLB'))) md[,(pair[2]):=num(get(pair[1]))]
for(z in c('Age','PMI')) md[,(z):=num(get(z))]
# Authoritative RIN: original RNA sample workbook, matched by NPID.
# DCL Patient ID is e.g. CWOW_NA20121, while cohort sample IDs can have
# additional group suffixes. The NPID component is used only if unique.
rna <- as.data.table(readxl::read_excel(RNA_SAMPLES_XLSX, sheet='Sheet1'))
if(!all(c('DCL Patient ID','RIN') %in% names(rna))) stop('RNA workbook needs DCL Patient ID and RIN columns')
rna[,npid_key:=toupper(trimws(as.character(`DCL Patient ID`)))]
rna[,RIN_workbook:=num(RIN)]
rna <- rna[!is.na(npid_key)&nzchar(npid_key)]
rin_conflicts <- rna[,.(n_rows=.N,n_distinct_values=uniqueN(RIN_workbook[is.finite(RIN_workbook)]),values=paste(unique(RIN_workbook[is.finite(RIN_workbook)]),collapse=';')),by=npid_key]
write_tab(rin_conflicts[n_rows>1 | n_distinct_values>1],'RIN_workbook_duplicate_ID_diagnostics.tsv')
if(any(rin_conflicts$n_distinct_values>1)) stop('Conflicting RIN values for identical DCL Patient ID; review RIN_workbook_duplicate_ID_diagnostics.tsv')
rna <- rna[,.(RIN_workbook=if(any(is.finite(RIN_workbook))) RIN_workbook[which(is.finite(RIN_workbook))[1L]] else NA_real_),by=npid_key]
md[,npid_rin_key:=toupper(sub('^CWOW_','',sample_id))]
md[,npid_rin_key:=sub('[-_].*$','',npid_rin_key)]
rna[,npid_key:=sub('^CWOW_','',npid_key)]
# Workbook IDs may include group suffixes (e.g. CWOW_NA19060-LBD-ATS).
# Strip suffixes consistently on BOTH sides of the NPID join.
rna[,npid_key:=sub('[-_].*$','',npid_key)]
if(anyDuplicated(rna$npid_key)) stop('Duplicate normalized RNA workbook NPID keys')
md[,RIN_cohort_original:=if('Tapestation_DIN_or_RIN'%in%names(md)) num(`Tapestation_DIN_or_RIN`) else NA_real_]
md[,RIN:=rna$RIN_workbook[match(npid_rin_key,rna$npid_key)]]
message('Workbook RIN coverage: ',sum(is.finite(md$RIN)),' / ',nrow(md))
if(sum(is.finite(md$RIN)) < 0.8*nrow(md)) stop('Workbook RIN match rate below 80%; inspect normalized NPID keys before QC.')
write_tab(md[,.(sample_id,npid_rin_key,RIN,RIN_cohort_original,RIN_difference=RIN-RIN_cohort_original)],'RIN_original_workbook_vs_cohort_metadata.tsv')
if(anyNA(md$RIN)) warning('Missing workbook RIN for ',sum(is.na(md$RIN)),' cohort metadata records; see RIN_original_workbook_vs_cohort_metadata.tsv. No fallback to cohort RIN is used.')
md[,SMRTcell_ID:=factor(as.character(SMRTcell_ID))]
md[,plot_id:=sub('^CWOW_','',sample_id)]
# Read each sample's TWO SMRT cell records from the GrandCanyon JSON.
# RGID (e.g. SMRT3496_L005_ATGACG_A57572) is a read-group ID, not a
# distinct SMRT cell by itself. The L005/L006 components distinguish runs/
# lanes within a shared SMRT label. Preserve raw RGID, BAM path, library ID,
# PacBio run directory, and inferred labels separately for reproducibility.
json_files <- list.files(JSON_DIR, pattern='\\.json$', full.names=TRUE)
if(!length(json_files)) stop('No sample JSON files found in ',JSON_DIR)
parse_json <- function(f) {
  j <- tryCatch(jsonlite::fromJSON(f,simplifyVector=FALSE),error=function(e) stop('Cannot parse ',f,': ',conditionMessage(e)))
  rec <- j$dataFiles
  if(is.null(rec)||!length(rec)) return(data.table(json_file=f,json_sample=tools::file_path_sans_ext(basename(f)),rgid=NA_character_,rgsm=NA_character_,rglb=NA_character_,bam_path=NA_character_,pacbio_run=NA_character_,smrt_label=NA_character_,smrt_lane=NA_character_,status='no_dataFiles'))
  rbindlist(lapply(rec,function(x) {
    rgid <- as.character(x$rgid %||% NA_character_)
    bam <- as.character(x$path %||% NA_character_)
    # Run directory is the first path segment matching PacBio rNNNN_... .
    seg <- strsplit(bam,'/',fixed=TRUE)[[1]]
    run <- seg[grepl('^r[0-9]+_',seg)]
    smrt <- regmatches(rgid,regexpr('SMRT[0-9]+',rgid))
    lane <- regmatches(rgid,regexpr('L[0-9]+',rgid))
    data.table(json_file=f,json_sample=as.character(j$project %||% tools::file_path_sans_ext(basename(f))),rgid=rgid,rgsm=as.character(x$rgsm %||% NA_character_),rglb=as.character(x$rglb %||% NA_character_),bam_path=bam,pacbio_run=if(length(run))run[1] else NA_character_,smrt_label=if(length(smrt)&&nzchar(smrt))smrt else NA_character_,smrt_lane=if(length(lane)&&nzchar(lane))lane else NA_character_,status=if(is.na(rgid))'missing_rgid' else 'ok')
  }),fill=TRUE)
}
`%||%` <- function(a,b) if(is.null(a)||length(a)==0L) b else a
json_rg <- rbindlist(lapply(json_files,parse_json),fill=TRUE)
json_rg[,sample_id:=fifelse(!is.na(rgsm)&nzchar(rgsm),rgsm,json_sample)]
json_rg[,sample_id:=paste0('CWOW_',sub('^CWOW_','',sample_id))]
# Match first by complete rgsm sample ID; fallback to the JSON project/NPID
# ONLY when the NPID identifies exactly one cohort donor.
json_rg[,npid_key:=sub('^CWOW_','',json_sample)]
md[,npid_key:=sub('^CWOW_','',sample_id)]
md[,npid_key:=sub('[-_].*$','',npid_key)]
json_rg[,npid_key:=sub('[-_].*$','',npid_key)]
npid_unique <- md[,.(n_cohort_samples=.N,cohort_sample_id=sample_id[1L]),by=npid_key]
json_rg <- merge(json_rg,npid_unique,by='npid_key',all.x=TRUE,sort=FALSE)
json_rg[!sample_id %in% md$sample_id & n_cohort_samples==1L,
        sample_id:=cohort_sample_id]
json_rg[,c('cohort_sample_id','n_cohort_samples'):=NULL]
# The original read-group identifiers and JSON filenames are retained in TSVs.

# RGSM is the complete donor-group sample ID; project/filename is only NPID.
# Export a durable long-format record of EVERY original SMRT/read group.
setcolorder(json_rg,c('sample_id',setdiff(names(json_rg),'sample_id')))
write_tab(json_rg,'SMRTcell_RGID_by_NPID_long.tsv')
json_rg <- json_rg[sample_id %in% md$sample_id]
if(anyDuplicated(json_rg[,.(sample_id,rgid)])) warning('Duplicate sample/RGID records detected; see long TSV')
json_summary <- json_rg[,.(n_datafiles=.N,n_distinct_rgid=uniqueN(rgid[!is.na(rgid)]),rgid_1=if(.N>=1)rgid[1] else NA_character_,rgid_2=if(.N>=2)rgid[2] else NA_character_,all_rgids=paste(sort(unique(na.omit(rgid))),collapse=';'),all_smrt_labels=paste(sort(unique(na.omit(smrt_label))),collapse=';'),all_smrt_lanes=paste(sort(unique(na.omit(smrt_lane))),collapse=';'),all_pacbio_runs=paste(sort(unique(na.omit(pacbio_run))),collapse=';'),all_libraries=paste(sort(unique(na.omit(rglb))),collapse=';')),by=sample_id]
md <- merge(md,json_summary,by='sample_id',all.x=TRUE,sort=FALSE)
if(anyNA(md$n_datafiles)) {
  missing_json <- md[is.na(n_datafiles),.(sample_id,NPID=if('NPID'%in%names(md)) as.character(NPID) else sub('^CWOW_','',sample_id),group_code)]
  write_tab(missing_json,'SMRTcell_missing_JSON_samples.tsv')
  warning('JSON missing or unmatched for ',nrow(missing_json),' count-matrix samples; see SMRTcell_missing_JSON_samples.tsv. This cannot be resolved without the correct JSON.')
} else {
  write_tab(data.table(sample_id=character(),NPID=character(),group_code=character()),'SMRTcell_missing_JSON_samples.tsv')
}
if(any(md$n_distinct_rgid!=2L,na.rm=TRUE)) warning('Some samples do not have exactly two distinct RGIDs; inspect SMRTcell_JSON_validation.tsv')
md[,json_validation:=fifelse(is.na(n_datafiles),'missing_json',fifelse(n_datafiles!=2L|n_distinct_rgid!=2L,'not_two_distinct_rgids','two_distinct_rgids'))]
write_tab(md[,.(sample_id,NPID=if('NPID'%in%names(md)) as.character(NPID) else sub('^CWOW_','',sample_id),group,sex,n_datafiles,n_distinct_rgid,rgid_1,rgid_2,all_rgids,all_smrt_labels,all_smrt_lanes,all_pacbio_runs,all_libraries,json_validation)],'SMRTcell_RGID_by_NPID_wide.tsv')
write_tab(md[,.(sample_id,json_validation,n_datafiles,n_distinct_rgid)],'SMRTcell_JSON_validation.tsv')
# Multiple cells per donor: record multi-membership incidence; this is NOT
# a conventional one-batch-per-sample factor for DESeq2.
incidence <- unique(json_rg[!is.na(smrt_label),.(sample_id,smrt_label)])
write_tab(incidence,'SMRTcell_sample_membership.tsv')
write_tab(incidence[,.(n_samples=uniqueN(sample_id)),by=smrt_label][order(-n_samples)],'SMRTcell_label_sample_counts.tsv')
# A one-to-one run-pair label is valid as a sample-level categorical descriptor,
# but can be high-cardinality or confounded; do not automatically fit it in DTE.
md[,run_pair:=factor(fifelse(is.na(all_pacbio_runs),NA_character_,all_pacbio_runs))]

# NPID_group ID: use NPID + group when NPID is available; otherwise sample ID.
if('NPID'%in%names(md)) {
  md[,NPID_group:=paste0(trimws(as.character(NPID)),'_',as.character(group))]
  md[is.na(NPID)|!nzchar(trimws(as.character(NPID))),NPID_group:=plot_id]
} else md[,NPID_group:=plot_id]
metadata_before_count_matching <- copy(md)
message('Reading moderate-confidence FL count matrix...')
raw <- fread(COUNTS,check.names=FALSE,showProgress=TRUE)
if(names(raw)[1]!='transcript_id') stop('First count column must be transcript_id')
ids <- as.character(raw$transcript_id)
if(anyDuplicated(ids)) stop('Duplicate transcript IDs')
samples <- names(raw)[-1L]
if(length(samples)!=274L) warning('Expected 274 samples; found ',length(samples))
md <- md[match(samples,sample_id)]
if(anyNA(md$sample_id)) stop('Count samples lack metadata: ',paste(head(samples[is.na(md$sample_id)],10),collapse=', '))
stopifnot(identical(md$sample_id,samples))
counts <- as.matrix(raw[,..samples]); rownames(counts) <- ids; rm(raw); gc()
if(anyNA(counts)||any(counts<0)||any(counts!=round(counts))) stop('Invalid FL count values')
storage.mode(counts) <- 'integer'
md[,`:=`(total_FL=colSums(counts),detected_transcripts=colSums(counts>0),detected_FL5=colSums(counts>=5))]
# Comprehensive reusable cohort metadata: original n276 CSV columns unchanged,
# plus authoritative workbook RIN, BOTH SMRT RGIDs, read-group/run/library
# descriptors, matching diagnostics, and expression QC. One row per original
# cohort record (including records not present in the 274-sample count matrix).
# Avoid overwriting original fields such as group, Gender, SMRTcell_ID, or
# Tapestation_DIN_or_RIN. The corrected RIN is an ADDITIONAL column.
aug_cols <- c('sample_id','npid_rin_key','RIN','RIN_cohort_original',
              'sex','group_code','NPID_group','n_datafiles','n_distinct_rgid',
              'rgid_1','rgid_2','all_rgids','all_smrt_labels','all_smrt_lanes',
              'all_pacbio_runs','all_libraries','json_validation',
              'total_FL','detected_transcripts','detected_FL5')
aug_cols <- intersect(aug_cols,names(md))
aug <- copy(md[,..aug_cols])
# Rebuild original cohort sample identifiers for a lossless one-to-one join.
md_original[,sample_id:=paste0('CWOW_',sub('^CWOW_','',trimws(as.character(`Study-Subject-ID`))))]
if(anyDuplicated(md_original$sample_id)) stop('Duplicate original cohort sample IDs')
# Keep original row order, and make count-matrix inclusion explicit.
aug[,in_count_matrix:=TRUE]
aug <- aug[match(md_original$sample_id,aug$sample_id)]
aug[,in_count_matrix:=fifelse(is.na(in_count_matrix),FALSE,in_count_matrix)]
# For original cohort rows outside the count matrix, retrieve workbook RIN
# and JSON-derived sequencing descriptors from pre-count matched metadata.
# This keeps the output faithful to ALL original n276 rows.
if(exists('metadata_before_count_matching')) {
  pre <- metadata_before_count_matching
  fill_cols <- intersect(c('RIN','RIN_cohort_original','npid_rin_key','sex',
                           'group_code','NPID_group','n_datafiles','n_distinct_rgid',
                           'rgid_1','rgid_2','all_rgids','all_smrt_labels',
                           'all_smrt_lanes','all_pacbio_runs','all_libraries',
                           'json_validation'),names(pre))
  missing_rows <- which(!aug$in_count_matrix)
  if(length(missing_rows)) {
    pre_idx <- match(md_original$sample_id[missing_rows],pre$sample_id)
    for(z in fill_cols) data.table::set(aug,i=missing_rows,j=z,value=pre[[z]][pre_idx])
  }
}
# NPID is retained from the original CSV. An explicit alias helps future joins.
aug[,RIN_source:='CWOW_Original_RNA_Samples.xlsx: RIN (DCL Patient ID)']
aug[,SMRT_source:='GrandCanyon JSON: dataFiles.rgid']
aug[,RIN_workbook:=RIN]
# Add new fields only; do not replace any source CSV columns.
new_cols <- setdiff(names(aug),names(md_original))
metadata_full <- cbind(copy(md_original),aug[,..new_cols])
stopifnot(nrow(metadata_full)==nrow(md_original),
          identical(metadata_full$sample_id,md_original$sample_id))
write_tab(metadata_full,'CWOW_n276_metadata_RIN_SMRTcell_RGID_augmented.tsv')
write_tab(metadata_full[in_count_matrix==TRUE],
          'CWOW_n274_DTE_metadata_RIN_SMRTcell_RGID_augmented.tsv')
write_tab(md,'matched_sample_metadata_and_yield.tsv')
vars <- c('RIN','Brain_weight','PMI','Age','Braak','Thal','CingLB','total_FL','detected_transcripts')
write_tab(data.table(variable=vars,n_available=vapply(md[,..vars],function(x)sum(is.finite(x)),integer(1)),n_missing=vapply(md[,..vars],function(x)sum(!is.finite(x)),integer(1))),'metadata_missingness.tsv')
write_tab(md[,.(n=.N),by=.(group,sex)],'group_by_sex.tsv')
write_tab(md[,.(n=.N),by=.(group,SMRTcell_ID)],'SMRTcell_by_group_confounding.tsv')
write_tab(unique(merge(incidence,md[,.(sample_id,group,sex)],by='sample_id'))[,.(n_samples=uniqueN(sample_id)),by=.(smrt_label,group)],'SMRTcell_label_by_group_confounding.tsv')
write_tab(md[,.(n=.N),by=.(run_pair,group)],'PacBio_run_pair_by_group_confounding.tsv')
# KW tests across five disease groups, overall and independently by sex.
# Every figure includes its actual p-value; all tests are consolidated in one TSV.
KW_RESULTS <- list()
kw_test <- function(z,y,stratum) {
  z <- z[is.finite(get(y)) & !is.na(group)]
  sizes <- z[,.(n=.N),by=group]
  if(nrow(z)<5L || uniqueN(z$group)<2L) return(data.table(variable=y,stratum=stratum,n=nrow(z),n_groups=uniqueN(z$group),H=NA_real_,df=NA_real_,p=NA_real_,group_n=paste(paste(sizes$group,sizes$n,sep=':'),collapse=';')))
  tst <- kruskal.test(z[[y]]~droplevels(z$group))
  data.table(variable=y,stratum=stratum,n=nrow(z),n_groups=uniqueN(z$group),H=unname(tst$statistic),df=unname(tst$parameter),p=tst$p.value,group_n=paste(paste(sizes$group,sizes$n,sep=':'),collapse=';'))
}
p_text <- function(p) ifelse(is.na(p),'KW p = NA',paste0('KW p = ',format.pval(p,digits=3,eps=1e-16)))
kw_violin <- function(d,y,ylab,stem) {
  z <- d[is.finite(get(y)) & !is.na(group)]
  if(nrow(z)<10L || uniqueN(z$group)<2L) return(invisible(NULL))
  overall <- kw_test(z,y,'All')
  female <- kw_test(z[sex=='Female'],y,'Female')
  male <- kw_test(z[sex=='Male'],y,'Male')
  KW_RESULTS[[length(KW_RESULTS)+1L]] <<- rbindlist(list(overall,female,male),fill=TRUE)
  p <- ggplot(z,aes(group,.data[[y]],fill=group))+
    geom_violin(alpha=.52,trim=FALSE,na.rm=TRUE)+
    geom_boxplot(width=.15,outlier.shape=NA,fill='white',alpha=.65)+
    geom_jitter(width=.11,size=1,alpha=.48)+scale_group_fill()+base_theme+
    theme(legend.position='none',plot.subtitle=element_text(face='bold',size=12))+
    labs(title=ylab,subtitle=p_text(overall$p),x=NULL,y=ylab)
  save_plot(p,stem,10,6)
  z2 <- z[!is.na(sex)]
  sex_labels <- rbindlist(list(female,male))[,.(sex=factor(stratum,levels=c('Female','Male')),label=p_text(p))]
  sex_labels[,`:=`(group=factor('Control',levels=GROUP_LABELS),y_position=Inf)]
  p2 <- ggplot(z2,aes(group,.data[[y]],fill=group))+
    geom_violin(alpha=.52,trim=FALSE,na.rm=TRUE)+
    geom_boxplot(width=.15,outlier.shape=NA,fill='white',alpha=.65)+
    geom_jitter(width=.12,size=.8,alpha=.45)+facet_wrap(~sex)+scale_group_fill()+
    geom_text(data=sex_labels,aes(x=group,y=y_position,label=label),inherit.aes=FALSE,
              hjust=0,vjust=1.4,fontface='bold',size=3.6)+
    scale_y_continuous(expand=expansion(mult=c(.05,.15)))+
    base_theme+theme(legend.position='none')+
    labs(title=paste(ylab,'by sex'),x=NULL,y=ylab)
  save_plot(p2,paste0(stem,'_by_sex'),13,6)
}
kw_violin(md,'RIN','RNA integrity (RIN)','Figure_01_RIN_by_group')
for(z in c('Brain_weight','PMI','Age','total_FL','detected_transcripts','Braak','Thal','CingLB'))
  kw_violin(md,z,c(Brain_weight='Brain weight (g)',PMI='Postmortem interval',Age='Age (years)',total_FL='Total FL support',detected_transcripts='Detected transcripts',Braak='Braak NFT stage',Thal='Thal amyloid phase',CingLB='Cingulate Lewy body burden')[[z]],paste0('Figure_07_KW_',z))
kw_all <- rbindlist(KW_RESULTS,fill=TRUE)
kw_all[,FDR_BH:=p.adjust(p,method='BH')]
write_tab(kw_all,'Kruskal_Wallis_all_variables_overall_and_by_sex.tsv')
for(z in c('total_FL','detected_transcripts','Brain_weight','PMI')) {
  d <- md[is.finite(RIN)&is.finite(get(z))]
  p <- ggplot(d,aes(RIN,.data[[z]],color=group))+geom_point(alpha=.7)+geom_smooth(method='lm',se=FALSE,linewidth=.6)+scale_group_color()+base_theme+labs(title=paste('RIN versus',z),x='RIN',y=if(z=='Brain_weight')'Brain weight (g)' else z)
  save_plot(p,paste0('Figure_02_RIN_vs_',z))
}
rin_cor <- rbindlist(lapply(c('total_FL','detected_transcripts','Brain_weight','PMI','Braak','Thal','CingLB'),function(z){d<-md[is.finite(RIN)&is.finite(get(z))];if(nrow(d)<5)return(data.table(variable=z,n=nrow(d),rho=NA_real_,p=NA_real_));v<-suppressWarnings(cor.test(d$RIN,d[[z]],method='spearman',exact=FALSE));data.table(variable=z,n=nrow(d),rho=unname(v$estimate),p=v$p.value)}))
rin_cor[,padj:=p.adjust(p,'BH')];write_tab(rin_cor,'RIN_spearman_associations.tsv')
# Global prevalence filter is for exploratory QC ONLY, not the DTE pairwise filter.
keep <- rowSums(counts>=MIN_FL)>=MIN_SAMPLES
if(sum(keep)<100L) stop('Too few transcripts retained for QC')
qc <- counts[keep,,drop=FALSE];rm(counts);gc()
message('QC transcripts: ',nrow(qc))
dds <- DESeqDataSetFromMatrix(qc,as.data.frame(md),design=~1)
ratio <- tryCatch(sizeFactors(estimateSizeFactors(dds,type='ratio')),error=function(e){message('Ratio size factors unavailable: ',conditionMessage(e));NULL})
pos <- sizeFactors(estimateSizeFactors(dds,type='poscounts'))
size_dt <- data.table(sample_id=samples,NPID_group=md$NPID_group,group=md$group,sex=md$sex,RIN=md$RIN,library_FL_QC=colSums(qc),size_factor_poscounts=pos)
if(!is.null(ratio)) size_dt[,size_factor_ratio:=ratio]
write_tab(size_dt,'normalization_size_factors.tsv')
# Robust outlier labeling in figures 03/04: points with extreme robust distance.
label_outliers <- function(d,columns,stem) {
  z <- copy(d)
  for(k in columns){v<-log10(z[[k]]);m<-median(v,na.rm=TRUE);s<-mad(v,constant=1.4826,na.rm=TRUE);z[,paste0('Z_',k):=if(is.finite(s)&&s>0) abs(v-m)/s else 0]}
  zs <- paste0('Z_',columns);z[,outlier_score:=do.call(pmax,c(.SD,list(na.rm=TRUE))),.SDcols=zs]
  z[,is_outlier:=outlier_score>=3.5]
  if(!any(z$is_outlier,na.rm=TRUE)) z[which.max(outlier_score),is_outlier:=TRUE]
  write_tab(z[is_outlier==TRUE,.(sample_id,NPID_group,outlier_score)],paste0(stem,'_outliers.tsv'))
  z
}
plot_norm <- function(d,x,y,stem,title,xlab,ylab) {
  z <- label_outliers(d,c(x,y),stem)
  p <- ggplot(z,aes(x=.data[[x]],y=.data[[y]],color=group))+geom_point(size=2,alpha=.85)+scale_x_log10()+scale_y_log10()+scale_group_color()+base_theme+labs(title=title,x=xlab,y=ylab)
  labs <- z[is_outlier==TRUE]
  if(nrow(labs)) {
    if(requireNamespace('ggrepel',quietly=TRUE)) p<-p+ggrepel::geom_text_repel(data=labs,aes(label=NPID_group),size=3,max.overlaps=Inf,show.legend=FALSE)
    else p<-p+geom_text(data=labs,aes(label=NPID_group),hjust=-.1,vjust=-.6,size=3,show.legend=FALSE)
  }
  save_plot(p,stem,10,7)
}
plot_norm(size_dt,'library_FL_QC','size_factor_poscounts','Figure_03_size_factor_vs_FL','Positive-count normalization versus FL library size','QC FL support','DESeq2 poscounts size factor')
if(!is.null(ratio)) plot_norm(size_dt,'size_factor_ratio','size_factor_poscounts','Figure_04_ratio_vs_poscounts','DESeq2 normalization comparison','Ratio size factor','Poscounts size factor')
# PCA on the most variable log2(poscounts-normalized FL+1) transcripts.
logmat <- log2(sweep(qc,2,pos,'/')+1)
rv <- rowVars(logmat);top <- order(rv,decreasing=TRUE)[seq_len(min(TOP_VARIABLE,length(rv)))]
mat <- t(logmat[top,,drop=FALSE]);mat <- scale(mat,center=TRUE,scale=FALSE)
pcs <- prcomp(mat,center=FALSE,scale.=FALSE,rank.=min(20L,nrow(mat)-1L,ncol(mat)))
explained <- 100*pcs$sdev^2/sum(pcs$sdev^2)
write_tab(data.table(PC=paste0('PC',seq_along(explained)),percent_variance=explained),'PCA_variance_explained.tsv')
score <- copy(md);score <- cbind(score,as.data.table(pcs$x))
write_tab(score,'sample_PCA_scores_and_metadata.tsv')
pca_plot <- function(z,stem) {
  d <- score[!is.na(PC1)&!is.na(PC2)]
  p <- ggplot(d,aes(PC1,PC2,color=.data[[z]]))+geom_point(size=2,alpha=.85)+base_theme+labs(title=paste('Expression PCA:',z),x=sprintf('PC1 (%.1f%%)',explained[1]),y=sprintf('PC2 (%.1f%%)',explained[2]),color=z)
  if(z=='group')p<-p+scale_group_color()
  if(z=='sex')p<-p+scale_color_manual(values=SEX_COLORS,drop=FALSE)
  if(z%in%c('RIN','Brain_weight','Braak','Thal','CingLB','Age','PMI'))p<-p+scale_color_viridis_c(na.value='grey80')
  save_plot(p,stem)
}
for(z in c('group','sex','SMRTcell_ID','run_pair','RIN','Brain_weight','Braak','Thal','CingLB','Age','PMI'))pca_plot(z,paste0('Figure_05_PCA_',z))
pc_associations <- rbindlist(lapply(paste0('PC',1:10),function(axis){rbindlist(lapply(vars,function(z){d<-score[is.finite(get(axis))&is.finite(get(z))];if(nrow(d)<8)return(data.table(PC=axis,variable=z,n=nrow(d),rho=NA_real_,p=NA_real_));v<-suppressWarnings(cor.test(d[[axis]],d[[z]],method='spearman',exact=FALSE));data.table(PC=axis,variable=z,n=nrow(d),rho=unname(v$estimate),p=v$p.value)}))}))
pc_associations[,padj:=p.adjust(p,'BH')];write_tab(pc_associations,'PCA_PC1_10_metadata_associations.tsv')
# Triangular clinical/neuropathology correlation heatmap matching the
# supplied example. This is a Spearman correlation plot, not canonical CCA.
# Disease group and sex are encoded as categorical integers only for this
# exploratory correlation visualization; correlations involving these codes
# are coding-dependent and must not be interpreted as invariant effect sizes.
clinical_cols <- c('Age','group','Braak','Thal','RIN','CingLB','Brain_weight','sex')
clinical_names <- c('Age','Disease group','Braak NFT','Thal amyloid','RIN','Cingulate LB','Brain weight (g)','Sex')
clinical <- copy(score[,..clinical_cols])
clinical[,group:=as.numeric(group)]
clinical[,sex:=as.numeric(sex)]
setnames(clinical,clinical_cols,clinical_names)
cm <- suppressWarnings(cor(as.matrix(clinical),use='pairwise.complete.obs',method='spearman'))
write_tab(data.table(variable=rownames(cm),as.data.table(cm)),'Clinical_neuropathology_Spearman_correlation_matrix.tsv')
heat <- as.data.table(as.table(cm));setnames(heat,c('row','column','rho'))
heat[,`:=`(row_i=match(as.character(row),clinical_names),col_i=match(as.character(column),clinical_names))]
heat <- heat[row_i<col_i]
heat[,`:=`(row=factor(as.character(row),levels=rev(clinical_names)),column=factor(as.character(column),levels=clinical_names))]
p <- ggplot(heat,aes(column,row,fill=rho))+
  geom_tile(color='white',linewidth=.7)+geom_text(aes(label=ifelse(is.na(rho),'NA',sprintf('%.2f',rho))),size=3.4)+
  scale_fill_gradient2(low='#377EB8',mid='white',high='#B2182B',midpoint=0,limits=c(-1,1),na.value='grey90')+
  coord_fixed()+theme_minimal(base_size=12)+
  theme(panel.grid=element_blank(),axis.title=element_blank(),axis.text.x=element_text(angle=50,hjust=0,vjust=0),axis.text.y=element_text(color='black'),legend.title=element_blank())+
  labs(title='Clinical and neuropathology correlations',subtitle='Pairwise-complete Spearman rho; categorical codes are exploratory')
save_plot(p,'Figure_10_clinical_neuropathology_triangular_correlations',11,8)
# CCA: PC1-10 versus clinical/pathology measurements; complete-case only.
cca_pc <- paste0('PC',1:10);cca_vars <- c('RIN','Brain_weight','PMI','Age','Braak','Thal','CingLB')
cc <- score[complete.cases(score[,c(cca_pc,cca_vars),with=FALSE])]
write_tab(data.table(n_total=nrow(score),n_complete_CCA=nrow(cc),variables=paste(cca_vars,collapse=';')),'CCA_complete_case_summary.tsv')
if(nrow(cc)>max(30L,2L*(length(cca_pc)+length(cca_vars)))) {
  X<-scale(as.matrix(cc[,..cca_pc]));Y<-scale(as.matrix(cc[,..cca_vars]))
  good_x<-apply(X,2,function(z)all(is.finite(z)));good_y<-apply(Y,2,function(z)all(is.finite(z)))
  X<-X[,good_x,drop=FALSE];Y<-Y[,good_y,drop=FALSE]
  if(ncol(X)>=2&&ncol(Y)>=2&&qr(X)$rank==ncol(X)&&qr(Y)$rank==ncol(Y)) {
    fit<-cancor(X,Y)
    U<-X%*%fit$xcoef[,1:2,drop=FALSE];V<-Y%*%fit$ycoef[,1:2,drop=FALSE]
    cc[,`:=`(CCA_expression_1=U[,1],CCA_expression_2=U[,2],CCA_metadata_1=V[,1],CCA_metadata_2=V[,2])]
    write_tab(cc,'CCA_sample_scores.tsv')
    write_tab(data.table(axis=seq_along(fit$cor),canonical_correlation=fit$cor),'CCA_correlations.tsv')
    write_tab(data.table(variable=colnames(Y),axis1=fit$ycoef[,1],axis2=fit$ycoef[,2]),'CCA_metadata_coefficients.tsv')
    p<-ggplot(cc,aes(CCA_expression_1,CCA_expression_2,color=group))+geom_point(size=2,alpha=.8)+scale_group_color()+base_theme+labs(title='CCA: expression PC1-10 versus clinical and neuropathology variables',subtitle=paste('Complete cases:',nrow(cc),'| exploratory, unadjusted'),x='Expression canonical variate 1',y='Expression canonical variate 2')
    save_plot(p,'Figure_06_CCA_expression_scores',10,7)
    # Pairwise PC-clinical correlations, annotated in a readable rectangular heatmap.
    cor_mat<-cor(X,Y,method='spearman');heat<-as.data.table(as.table(cor_mat));setnames(heat,c('PC','variable','rho'))
    p<-ggplot(heat,aes(variable,PC,fill=rho))+geom_tile(color='white')+geom_text(aes(label=sprintf('%.2f',rho)),size=3)+scale_fill_gradient2(low='#3D75AD',mid='white',high='#C35D57',midpoint=0,limits=c(-1,1))+base_theme+labs(title='Spearman correlations: expression PC1-10 and clinical measures',x=NULL,y=NULL)
    save_plot(p,'Figure_08_PC1_10_clinical_correlation_heatmap',11,7)
  }else warning('CCA skipped: rank-deficient input')
}else warning('CCA skipped: insufficient complete cases')
# Per-transcript variance partition: fixed effects (variancePartition::fitExtractVarPartModel).
# Restrict to highly variable transcripts for manageable computation and output.
if(RUN_VARPART) {
  if(!requireNamespace('variancePartition',quietly=TRUE)) {
    warning('Install variancePartition to generate variance partition plots; skipping')
  } else {
    vp_vars<-c('group','sex','RIN','Age','PMI','Brain_weight')
    vp_md<-as.data.frame(md)
    good<-complete.cases(vp_md[,vp_vars,drop=FALSE])
    vmd<-droplevels(vp_md[good,,drop=FALSE])
    # Enforce sample identifiers on both sides of variancePartition.
    if(anyDuplicated(vmd$sample_id)) stop('Duplicate variance-partition sample IDs')
    rownames(vmd) <- vmd$sample_id
    vm<-logmat[order(rv,decreasing=TRUE)[seq_len(min(VARPART_TRANSCRIPTS,length(rv)))],good,drop=FALSE]
    if(!identical(colnames(vm),rownames(vmd))) {
      if(!setequal(colnames(vm),rownames(vmd))) stop('Variance-partition sample identifiers do not match metadata')
      vmd <- vmd[colnames(vm),,drop=FALSE]
    }
    stopifnot(identical(colnames(vm),rownames(vmd)))
    write_tab(data.table(expression_sample=colnames(vm),metadata_sample=rownames(vmd),matched=colnames(vm)==rownames(vmd)),
              'variance_partition_sample_alignment.tsv')
    if(ncol(vm)>30&&nlevels(vmd$group)>1&&nlevels(vmd$sex)>1) {
      # Remove invariant numeric predictors to prevent singular model fitting.
      numeric_vars<-c('RIN','Age','PMI','Brain_weight')
      active<-numeric_vars[vapply(vmd[,numeric_vars,drop=FALSE],sd,numeric(1),na.rm=TRUE)>0]
      form<-as.formula(paste('~',paste(c('group','sex',active),collapse=' + ')))
      message('Variance partition model: ',deparse(form))
      vp<-tryCatch(variancePartition::fitExtractVarPartModel(vm,form,vmd),error=function(e){warning('Variance partition failed: ',conditionMessage(e));NULL})
      if(!is.null(vp)) {
        write_tab(data.table(transcript_id=rownames(vp),as.data.table(vp)),'variance_partition_transcripts.tsv')
        vp_long<-melt(data.table(transcript_id=rownames(vp),as.data.table(vp)),id.vars='transcript_id',variable.name='term',value.name='fraction')
        p<-ggplot(vp_long,aes(term,100*fraction,fill=term))+geom_violin(trim=FALSE,alpha=.65)+geom_boxplot(width=.13,fill='white',outlier.shape=NA)+base_theme+theme(legend.position='none')+labs(title='Per-transcript variance partition',subtitle=paste(nrow(vp),'high-variance transcripts; complete-case fixed-effects model'),x=NULL,y='Variance explained (%)')
        save_plot(p,'Figure_09_variance_partition',11,7)
        write_tab(vp_long[,.(median_percent=100*median(fraction,na.rm=TRUE),mean_percent=100*mean(fraction,na.rm=TRUE)),by=term],'variance_partition_summary.tsv')
      }
    } else warning('Variance partition skipped: insufficient samples or group/sex levels')
  }
}
writeLines(c(paste('Count matrix:',COUNTS),paste('Metadata:',METADATA),paste('RIN workbook:',RNA_SAMPLES_XLSX),paste('JSON directory:',JSON_DIR),paste('Output:',OUT),'SMRT read groups: 2 per donor expected; raw RGIDs and PacBio runs exported as long and wide TSV; multi-membership is not a simple DESeq2 batch factor','PCA: log2(poscounts-normalized FL+1), top-variable transcripts','RIN is exclusively from the original RNA samples workbook (DCL Patient ID); missing values remain NA','KW overall and sex-stratified tests with BH-adjusted p-values in consolidated TSV','Triangular clinical Spearman correlation heatmap is distinct from canonical CCA','CCA exploratory; pathology not adjusted from DTE',paste('Variance partition fixed effects, top',VARPART_TRANSCRIPTS,'high-variance transcripts'), 'Augmented metadata preserves all original n276 columns and adds workbook RIN plus JSON RGIDs'),file.path(OUT,'analysis_notes.txt'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
message('Complete: ',OUT)
