# =============================================================================
# CWOW | SMRT sequencing-run confounding and expression-PC associations
# Position: after 03_PreDTE_covariate_assessment_RIN_JOIN_FIXED.R;
#           BEFORE selecting covariates for pairwise DESeq2 DTE.
# Directly source in Posit Workbench. Read-only inputs; no DTE rerun.
# INPUTS: augmented n274 metadata, sample PCA scores, SMRT sample membership.
# OUTPUTS: TSV diagnostics, PDF/PNG figures, interpretation checklist.
# IMPORTANT: 2 RGIDs per biological donor are NOT independent samples.
#           'run_pair' is a sample-level composite technical descriptor;
#           it is NOT proof of an independent sequencing batch effect.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
ROOT <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW'
QC_DIR <- file.path(ROOT,'results/QC_technical_RIN_normalization_CCA')
METADATA <- file.path(QC_DIR,'CWOW_n274_DTE_metadata_RIN_SMRTcell_RGID_augmented.tsv')
PCA_FILE <- file.path(QC_DIR,'sample_PCA_scores_and_metadata.tsv')
MEMBERSHIP_FILE <- file.path(QC_DIR,'SMRTcell_sample_membership.tsv')
OUT <- file.path(QC_DIR,'SMRT_run_confounding_PC_assessment')
N_PCS <- 10L
MIN_SAMPLES_PER_RUN_PAIR <- 3L # reporting flag; do not discard small batches
GROUP_COLORS <- c('Control'='#377EB8','AD(AT)'='#E41A1C','LBD(S)'='#6A3D9A',
                  'LBD(AS)'='#9E9AC8','LBD(ATS)'='#636363')
for(f in c(METADATA,PCA_FILE,MEMBERSHIP_FILE)) if(!file.exists(f)) stop('Missing input: ',f)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write_tab <- function(x,name) fwrite(as.data.table(x),file.path(OUT,name),sep='\t',na='NA')
save_plot <- function(p,name,w=10,h=7) {
  ggsave(file.path(OUT,paste0(name,'.png')),p,width=w,height=h,dpi=300,bg='white',limitsize=FALSE)
  ggsave(file.path(OUT,paste0(name,'.pdf')),p,width=w,height=h,bg='white',limitsize=FALSE)
  message('Saved ',name)
}
fmt <- function(x) ifelse(is.na(x),'NA',formatC(x,format='fg',digits=3))
message('Loading existing sample-level QC metadata and expression PCs...')
md <- fread(METADATA,check.names=FALSE)
pc <- fread(PCA_FILE,check.names=FALSE)
member <- fread(MEMBERSHIP_FILE,check.names=FALSE)
needed <- c('sample_id','group','sex','Age','RIN','rgid_1','rgid_2','all_smrt_labels','all_pacbio_runs')
if(!all(needed %in% names(md))) stop('Missing metadata columns: ',paste(setdiff(needed,names(md)),collapse=', '))
if(!all(c('sample_id','smrt_label') %in% names(member))) stop('Membership table needs sample_id and smrt_label')
if(!'sample_id' %in% names(pc)) stop('PCA table needs sample_id')
if(anyDuplicated(md$sample_id)||anyDuplicated(pc$sample_id)) stop('Duplicate sample_id in metadata/PCA')
pcs <- paste0('PC',seq_len(N_PCS))
if(!all(pcs %in% names(pc))) stop('PCA scores missing: ',paste(setdiff(pcs,names(pc)),collapse=', '))
# Use authoritative augmented metadata fields, not similarly named PCA-file fields.
score <- merge(md,pc[,c('sample_id',pcs),with=FALSE],by='sample_id',all.x=TRUE,sort=FALSE)
if(anyNA(score$PC1)) stop('PCA scores absent for some samples; rerun the QC PCA on matching n274 samples')
score[,group:=factor(as.character(group),levels=c('CONTROL','AD','LBD_S','LBD_AS','LBD_ATS'),
                     labels=names(GROUP_COLORS))]
score[,sex:=factor(as.character(sex))]
score[,Age:=suppressWarnings(as.numeric(Age))]
score[,RIN:=suppressWarnings(as.numeric(RIN))]
score[,Age_centered:=Age-mean(Age,na.rm=TRUE)]
# Run pairs are sorted from the two read-group-derived run IDs already in the metadata.
clean_pair <- function(x) {
  if(is.na(x)||!nzchar(trimws(x))) return(NA_character_)
  paste(sort(unique(trimws(strsplit(x,';',fixed=TRUE)[[1]]))),collapse=';')
}
score[,run_pair:=vapply(as.character(all_pacbio_runs),clean_pair,character(1))]
score[,smrt_pair:=vapply(as.character(all_smrt_labels),clean_pair,character(1))]
score[,run_pair:=factor(run_pair)]
score[,smrt_pair:=factor(smrt_pair)]
write_tab(score[,c('sample_id','NPID','NPID_group','group','sex','Age','Age_centered','RIN',
                    'rgid_1','rgid_2','all_smrt_labels','all_pacbio_runs','run_pair',pcs),with=FALSE],
          '01_sample_run_pair_and_PCA.tsv')
# Validation: exactly two read groups, and run membership is not expanded to pseudo-replicates.
score[,n_rgids:=as.integer(!is.na(rgid_1)&nzchar(rgid_1)) + as.integer(!is.na(rgid_2)&nzchar(rgid_2))]
write_tab(score[,.(sample_id,NPID,group,rgid_1,rgid_2,n_rgids,run_pair,smrt_pair)],
          '02_run_membership_validation.tsv')
if(any(score$n_rgids!=2L)) warning('Some samples lack two RGIDs; inspect 02_run_membership_validation.tsv')
# Composite run-pair distribution across disease groups.
run_ct <- dcast(score[,.(n=.N),by=.(run_pair,group)],run_pair~group,value.var='n',fill=0)
write_tab(run_ct,'03_run_pair_by_group_counts.tsv')
run_n <- score[,.(n_samples=.N,n_groups=uniqueN(group),
                   groups=paste(sort(unique(as.character(group))),collapse=';')),
               by=run_pair][order(n_samples)]
run_n[,small_batch:=n_samples<MIN_SAMPLES_PER_RUN_PAIR]
write_tab(run_n,'04_run_pair_sample_sizes.tsv')
# Multi-membership: one sample can appear in multiple SMRT labels; use for descriptive QC only.
member <- unique(member[,.(sample_id=as.character(sample_id),smrt_label=as.character(smrt_label))])
member <- merge(member,score[,.(sample_id,group)],by='sample_id',all.x=TRUE)
write_tab(member[,.(n_samples=uniqueN(sample_id)),by=.(smrt_label,group)],
          '05_SMRT_label_by_group_membership_long.tsv')
# Chi-squared is descriptive; sparse contingency cells make asymptotic p unreliable.
ct <- table(score$group,score$run_pair)
chisq <- suppressWarnings(chisq.test(ct,simulate.p.value=TRUE,B=9999))
cramer_v <- sqrt(unname(sum((ct-chisq$expected)^2/chisq$expected))/
                   (sum(ct)*min(nrow(ct)-1,ncol(ct)-1)))
write_tab(data.table(test='Group vs run_pair; Monte Carlo chi-square',n=sum(ct),
                     statistic=unname(chisq$statistic),p_value=chisq$p.value,
                     cramers_v=cramer_v,n_groups=nrow(ct),n_run_pairs=ncol(ct)),
          '06_group_run_pair_confounding_test.tsv')
write_tab(as.data.table(as.table(ct)),'07_group_run_pair_contingency_long.tsv')
plot_counts <- as.data.table(as.table(ct));setnames(plot_counts,c('group','run_pair','N'))
p <- ggplot(plot_counts,aes(x=run_pair,y=group,fill=N))+geom_tile(color='white')+
  geom_text(aes(label=ifelse(N>0,N,'')),size=3)+
  scale_fill_gradient(low='white',high='#377EB8')+theme_minimal(base_size=11)+
  theme(axis.text.x=element_text(angle=60,hjust=1),panel.grid=element_blank())+
  labs(title='Sequencing run-pair representation across neuropathology groups',
       subtitle='One observation per donor; two RGIDs retained in metadata',x='Run pair',y=NULL)
save_plot(p,'Figure_01_group_by_run_pair_heatmap',13,6)
p <- ggplot(score,aes(run_pair,fill=group))+geom_bar(position='fill')+
  scale_fill_manual(values=GROUP_COLORS,drop=FALSE)+scale_y_continuous(labels=function(x)paste0(round(x*100),'%'))+
  theme_classic(base_size=11)+theme(axis.text.x=element_text(angle=60,hjust=1))+
  labs(title='Disease-group composition within sequencing run pairs',x='Run pair',y='Percentage of donors',fill='Group')
save_plot(p,'Figure_02_group_composition_by_run_pair',13,6)
# Model estimability and confounding: compare covariate models on the SAME complete cases.
cc <- score[complete.cases(score[,.(group,sex,Age_centered,RIN,run_pair)])]
if(nrow(cc)<30) stop('Too few complete cases for age/sex/RIN/run-pair models')
cc[,run_pair:=droplevels(run_pair)];cc[,group:=droplevels(group)];cc[,sex:=droplevels(sex)]
forms <- list(
  baseline=~sex+Age_centered+RIN,
  disease=~sex+Age_centered+RIN+group,
  run_only=~sex+Age_centered+RIN+run_pair,
  disease_and_run=~sex+Age_centered+RIN+group+run_pair
)
model_diag <- rbindlist(lapply(names(forms),function(nm){
  mm <- model.matrix(forms[[nm]],as.data.frame(cc))
  q <- qr(mm)
  data.table(model=nm,n_samples=nrow(mm),n_columns=ncol(mm),rank=q$rank,
             residual_df=nrow(mm)-q$rank,full_rank=q$rank==ncol(mm),
             aliased_columns=if(q$rank<ncol(mm))paste(colnames(mm)[q$pivot[(q$rank+1):ncol(mm)]],collapse=';') else '',
             condition_number=tryCatch(kappa(mm),error=function(e)NA_real_))
}),fill=TRUE)
write_tab(model_diag,'08_DESeq2_design_rank_and_estimability.tsv')
# Is run_pair distinguishable from group? Incremental model rank tells us how many
# estimable dimensions each term contributes after accounting for the other.
rank_of <- function(f) qr(model.matrix(f,as.data.frame(cc)))$rank
increment <- data.table(n_complete=nrow(cc),
  group_df_after_run=rank_of(forms$disease_and_run)-rank_of(forms$run_only),
  run_df_after_group=rank_of(forms$disease_and_run)-rank_of(forms$disease),
  n_group_levels=nlevels(cc$group),n_run_pair_levels=nlevels(cc$run_pair),
  n_sex_levels=nlevels(cc$sex))
write_tab(increment,'09_incremental_rank_group_and_run.tsv')
if(!model_diag[model=='disease_and_run',full_rank]) {
  warning('Combined disease + run-pair design is rank deficient: DO NOT use it unchanged in DESeq2')
}
# Formal association of run pair with PC1-PC10 after controlling group, sex,
# centered age and RIN. Nested-model F tests are exploratory and unadjusted
# for transcript-level inference. No assumption that each RGID is a donor.
assoc <- rbindlist(lapply(pcs,function(axis){
  d <- cc[is.finite(get(axis))]
  if(nrow(d)<30)return(NULL)
  fit0 <- lm(reformulate(c('sex','Age_centered','RIN','group'),response=axis),data=d)
  fit1 <- lm(reformulate(c('sex','Age_centered','RIN','group','run_pair'),response=axis),data=d)
  fit2 <- lm(reformulate(c('sex','Age_centered','RIN','run_pair'),response=axis),data=d)
  an_run <- anova(fit0,fit1);an_group <- anova(fit2,fit1)
  r20 <- summary(fit0)$r.squared;r21 <- summary(fit1)$r.squared
  data.table(PC=axis,n=nrow(d),R2_without_run=r20,R2_with_run=r21,
             delta_R2_run=r21-r20,
             run_df=an_run$Df[2],run_F=an_run$F[2],run_p=an_run$`Pr(>F)`[2],
             group_df=an_group$Df[2],group_F=an_group$F[2],group_p=an_group$`Pr(>F)`[2])
}),fill=TRUE)
if(!nrow(assoc))stop('No expression PC association results')
assoc[,run_p_BH:=p.adjust(run_p,method='BH')]
assoc[,group_p_BH:=p.adjust(group_p,method='BH')]
write_tab(assoc,'10_PC1_10_run_and_group_partial_F_tests.tsv')
plot_a <- melt(assoc[,.(PC,delta_R2_run)],id.vars='PC',variable.name='metric',value.name='delta_R2')
p <- ggplot(plot_a,aes(x=PC,y=delta_R2,fill=PC))+geom_col(show.legend=FALSE)+
  geom_text(data=assoc,aes(x=PC,y=delta_R2_run,label=paste0('FDR=',fmt(run_p_BH))),
            inherit.aes=FALSE,vjust=-0.3,size=3)+
  coord_cartesian(ylim=c(0,max(assoc$delta_R2_run,na.rm=TRUE)*1.25+0.005))+
  theme_classic(base_size=11)+labs(title='Incremental variance in expression PCs associated with run pair',
    subtitle='Nested linear models adjusting for disease group, sex, centered age, and RIN',
    x='Expression PC',y='Change in R-squared when run pair is added')
save_plot(p,'Figure_03_PC1_10_incremental_run_R2',11,6)
# PCA plots use run pairs as a descriptive legend; groups use the fixed palette.
for(v in c('run_pair','group')) {
  p <- ggplot(score,aes(PC1,PC2,color=.data[[v]]))+geom_point(size=2,alpha=.85)+
    theme_classic(base_size=12)+labs(title=paste('Expression PCA by',v),
    x='PC1',y='PC2',color=if(v=='group')'Group' else 'Run pair')
  if(v=='group')p<-p+scale_color_manual(values=GROUP_COLORS,drop=FALSE)
  save_plot(p,paste0('Figure_04_PCA_',v),if(v=='run_pair')12 else 9,7)
}
# Quantify overlap of run-pair design with each disease group. No automatic
# batch selection: full-rank alone is not proof of a scientifically useful batch.
write_tab(cc[,.(n=.N,median_RIN=median(RIN,na.rm=TRUE),mean_age=mean(Age,na.rm=TRUE)),
             by=.(run_pair,group)][order(run_pair,group)],'11_run_pair_group_age_RIN.tsv')
notes <- c(
  paste('Total metadata donors:',nrow(score)),
  paste('Complete-case donors:',nrow(cc)),
  paste('Unique run-pair categories:',nlevels(cc$run_pair)),
  paste('Cramers V (group vs run pair):',fmt(cramer_v)),
  paste('Monte Carlo association p:',fmt(chisq$p.value)),
  paste('Combined DESeq2 design full rank:',model_diag[model=='disease_and_run',full_rank]),
  paste('Run dimensions after disease:',increment$run_df_after_group),
  paste('Disease dimensions after run:',increment$group_df_after_run),
  'INTERPRETATION: run-pair association with PCs does not prove a causal batch effect.',
  'Do not add two RGIDs as separate donor-level rows; each donor contributes one expression profile.',
  'Run pairs with few donors and groups represented in only one run require careful review.',
  'If group and run are perfectly confounded, no fixed-effects model can separate them.',
  'The primary DTE decision requires balancing estimability, effect size, power, and known experimental processing.'
)
writeLines(notes,file.path(OUT,'12_analysis_summary_and_cautions.txt'))
message('Complete. Results: ',OUT)
