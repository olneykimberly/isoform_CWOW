# ============================================================================
# CWOW FL audit: sample and disease-group plots
# Source directly in Posit Workbench. Does not modify audit or DTE inputs.
# Mapped FL = SQANTI3 full-length support successfully assigned to merged IDs.
# ============================================================================

# USER SETTINGS --------------------------------------------------------------
BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
AUDIT_DIR <- file.path(BASE, '05_DTE_DTU', 'QC_IsoMatch_FL_reconciliation_disk_backed')
INPUT <- file.path(AUDIT_DIR, 'sample_reconciliation.tsv')
OUT <- file.path(AUDIT_DIR, 'Sample_Reconciliation_Figures')
DPI <- 300
# ---------------------------------------------------------------------------

for (pkg in c('data.table', 'ggplot2', 'scales')) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop('Install missing package: ',pkg)
}
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
if (!file.exists(INPUT)) stop('Missing audit results: ',INPUT)

GROUP_LEVELS <- c('Control','AD(AT)','LBD(S)','LBD(AS)','LBD(ATS)')
GROUP_COLORS <- c('Control'='#377EB8','AD(AT)'='#E41A1C',
                  'LBD(S)'='#6A3D9A','LBD(AS)'='#9E9AC8','LBD(ATS)'='#636363')

log_file <- file.path(OUT, paste0('sample_reconciliation_plots_',format(Sys.time(),'%Y%m%d_%H%M%S'),'.log'))
log_con <- file(log_file,'wt')
log_msg <- function(...) {
  z <- paste0(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),' | ',paste0(...,collapse=''))
  cat(z,'\n'); cat(z,'\n',file=log_con); flush(log_con)
}

save_plot <- function(p, stem, width, height) {
  for (ext in c('png','pdf')) {
    f <- file.path(OUT,paste0(stem,'.',ext))
    ggplot2::ggsave(f,plot=p,width=width,height=height,units='in',
                    dpi=if(ext=='png') DPI else 300,limitsize=FALSE,bg='white')
    log_msg('Saved: ',f)
  }
}

tryCatch({
  d <- fread(INPUT)
  needed <- c('sample_id','source_total_FL','source_unmapped_FL',
              'mapped_expected_FL','observed_total_FL')
  if (!all(needed %in% names(d))) stop('Missing columns: ',paste(setdiff(needed,names(d)),collapse=', '))
  if (anyDuplicated(d$sample_id)) stop('Duplicate sample_id in reconciliation file')
  # Sample labels from the audit include disease suffixes, e.g. CWOW_NA00118-LBD-AS.
  # Use suffix anchoring so that LBD-AS is never confused with LBD-ATS.
  x <- toupper(trimws(d$sample_id))
  d[, group := fcase(
    grepl('[-_]CONTROL$',x), 'Control',
    grepl('[-_]AD(?:[-_]AT)?$',x), 'AD(AT)',
    grepl('[-_]LBD[-_]ATS$',x), 'LBD(ATS)',
    grepl('[-_]LBD[-_]AS$',x), 'LBD(AS)',
    grepl('[-_]LBD[-_]S$',x), 'LBD(S)',
    default=NA_character_)]
  if (anyNA(d$group)) stop('Cannot determine disease group from sample IDs: ',
                           paste(head(d[is.na(group),sample_id],15),collapse=', '))
  d[, group := factor(group,levels=GROUP_LEVELS)]
  d[, mapped_FL := as.numeric(mapped_expected_FL)]
  d[, unmapped_FL := as.numeric(source_unmapped_FL)]
  d[, mapping_pct := fifelse(source_total_FL>0,
                             100*mapped_FL/source_total_FL,NA_real_)]
  d[, reconciliation_pct := fifelse(mapped_FL>0,
                                   100*observed_total_FL/mapped_FL,NA_real_)]
  if (any(!is.finite(d$mapped_FL))) stop('Nonfinite mapped FL counts')
  setorder(d,group,-mapped_FL,sample_id)
  d[, sample_order := seq_len(.N)]
  d[, within_group_order := seq_len(.N),by=group]
  fwrite(d,file.path(OUT,'sample_reconciliation_with_groups.tsv'),sep='\t')
  summary <- d[,.(n_samples=.N,
                  total_mapped_FL=sum(mapped_FL),
                  median_mapped_FL=median(mapped_FL),
                  mean_mapped_FL=mean(mapped_FL),
                  min_mapped_FL=min(mapped_FL),
                  max_mapped_FL=max(mapped_FL),
                  median_mapping_pct=median(mapping_pct,na.rm=TRUE)),by=group]
  fwrite(summary,file.path(OUT,'mapped_FL_by_disease_group_summary.tsv'),sep='\t')
  log_msg('Read ',nrow(d),' samples; group counts: ',
          paste(paste(summary$group,summary$n_samples,sep='='),collapse=', '))

  fmt <- scales::label_comma(accuracy=1)
  base_theme <- theme_classic(base_size=12)+
    theme(plot.title=element_text(face='bold',size=14),
          plot.subtitle=element_text(size=10),
          legend.position='none')
  cols <- scale_fill_manual(values=GROUP_COLORS,drop=FALSE)

  # 1. All samples; ordered by disease group and descending mapped FL within group.
  p1 <- ggplot(d,aes(x=sample_order,y=mapped_FL,fill=group))+
    geom_col(width=.95)+cols+
    scale_y_continuous(labels=fmt,expand=expansion(mult=c(0,.05)))+
    scale_x_continuous(breaks=NULL,expand=expansion(mult=c(.002,.002)))+
    labs(title='Mapped full-length support by sample',
         subtitle='Samples ordered by neuropathological group and mapped FL within group',
         x='Samples (n = 274)',y='Mapped SQANTI3 FL counts')+base_theme+
    theme(axis.line.x=element_blank(),axis.ticks.x=element_blank())
  # Draw separators at group boundaries, with a group key below the plot.
  ends <- d[,.(end=max(sample_order)),by=group]$end
  p1 <- p1+geom_vline(xintercept=head(ends,-1)+.5,color='white',linewidth=.35)
  save_plot(p1,'Figure_01_Mapped_FL_All_Samples_Bar',16,6)

  # 2. Boxplots with sample dots: group differences and outliers visible.
  p2 <- ggplot(d,aes(x=group,y=mapped_FL,fill=group))+
    geom_boxplot(width=.63,outlier.shape=NA,alpha=.85,color='#303030',linewidth=.5)+
    geom_jitter(width=.17,height=0,size=1.35,alpha=.58,color='#252525')+
    cols+scale_y_continuous(labels=fmt,expand=expansion(mult=c(.03,.08)))+
    labs(title='Mapped full-length support by disease group',
         subtitle='Each point represents one sample; boxes show median and interquartile range',
         x=NULL,y='Mapped SQANTI3 FL counts')+base_theme+
    theme(axis.text.x=element_text(size=11))
  save_plot(p2,'Figure_02_Mapped_FL_Group_Boxplot',10,6)

  # 3. Sample bars faceted by group, free x-scale to accommodate different n.
  # Keep common y-scale to allow between-group comparison.
  p3 <- ggplot(d,aes(x=within_group_order,y=mapped_FL,fill=group))+
    geom_col(width=.92)+cols+
    facet_wrap(~group,ncol=1,scales='free_x')+
    scale_y_continuous(labels=fmt,expand=expansion(mult=c(0,.06)))+
    scale_x_continuous(breaks=NULL,expand=expansion(mult=c(.002,.002)))+
    labs(title='Mapped full-length support for each sample, by disease group',
         subtitle='Samples sorted from highest to lowest mapped FL within each group; shared y-axis scale',
         x='Individual samples (sorted within group)',y='Mapped SQANTI3 FL counts')+base_theme+
    theme(strip.background=element_rect(fill='#F2F2F2',color=NA),
          strip.text=element_text(face='bold',hjust=0),
          panel.spacing.y=grid::unit(.7,'lines'),
          axis.line.x=element_blank(),axis.ticks.x=element_blank())
  save_plot(p3,'Figure_03_Mapped_FL_Sample_Bars_Faceted',14,13)

  # 4. Diagnostic: percentage of source FL successfully mapped to merged IDs.
  p4 <- ggplot(d,aes(x=group,y=mapping_pct,fill=group))+
    geom_boxplot(width=.63,outlier.shape=NA,alpha=.85,color='#303030')+
    geom_jitter(width=.17,height=0,size=1.35,alpha=.58,color='#252525')+
    cols+scale_y_continuous(labels=function(z) paste0(z,'%'))+
    labs(title='Fraction of SQANTI3 FL support mapped to IsoMatch IDs',
         subtitle='Mapped expected FL / total source FL × 100; this is not the observed-versus-expected agreement',
         x=NULL,y='Source FL successfully mapped (%)')+base_theme
  save_plot(p4,'Figure_04_Mapping_Percentage_Group_Boxplot',10,6)

  # 5. Diagnostic: observed versus expected FL agreement, in counts.
  p5 <- ggplot(d,aes(x=mapped_FL,y=observed_total_FL,color=group))+
    geom_abline(slope=1,intercept=0,linetype='dashed',color='#777777')+
    geom_point(alpha=.8,size=2)+
    scale_color_manual(values=GROUP_COLORS,drop=FALSE)+
    scale_x_continuous(labels=fmt)+scale_y_continuous(labels=fmt)+
    coord_equal()+
    labs(title='FL count reconciliation across samples',
         subtitle='Dashed line indicates perfect agreement',
         x='Expected mapped SQANTI3 FL',y='Observed IsoMatch FL',color='Group')+
    base_theme+theme(legend.position='right')
  save_plot(p5,'Figure_05_Expected_vs_Observed_FL',8,7)

  log_msg('SUCCESS: all plots and summaries saved to ',OUT)
},error=function(e){log_msg('ERROR: ',conditionMessage(e));stop(e)},
finally={close(log_con)})
