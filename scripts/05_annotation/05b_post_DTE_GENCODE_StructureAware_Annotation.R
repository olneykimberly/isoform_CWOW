# CWOW: structure-aware GENCODE v38 annotation of IsoMatch transcripts
# Source directly in Posit Workbench. Does not rerun DTE or overwrite its results.
# Annotates full and significant DTE files; does not modify original results.
# A gene symbol is emitted only for a uniquely supported structural gene assignment.
# This is a conservative structural annotation, not SQANTI3 reclassification.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})


BASE <- '/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/results/isomatch_outputs'
SCENARIO <- 'moderate'
ISOMATCH_GTF <- file.path(BASE, '04_confidence_filter_rescuedLR_274', SCENARIO,
                           'CWOW_274_rescuedLR.confidence_filtered.gtf.gz')
GENCODE_GTF <- '/tgen_labs/jfryer/projects/references/human/GRCh38/gencode.v38.annotation.gtf'
DTE_DIR <- file.path(BASE, '05_DTE_DTU', 'DTE_DESeq2', SCENARIO)
OUT_DIR <- file.path(BASE, '05_DTE_DTU', 'GENCODE_structure_aware', SCENARIO)
RESULT_PATTERN <- '^DTE_(results|FDR05)\\.tsv(\\.gz)?$'  # significant results; see file discovery below
# Change RESULT_PATTERN to '^DTE_FDR05\\.tsv(\\.gz)?$' for significant results only.
MIN_EXON_OVERLAP <- 0.80       # fraction of IsoMatch exonic bases overlapping gene exons
MIN_WINNER_MARGIN <- 0.10     # fraction separating strongest and second-best candidates
CHUNK_LINES <- 100000L

stopifnot(file.exists(ISOMATCH_GTF), file.exists(GENCODE_GTF))
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

attr_get <- function(x, key) {
  pat <- paste0('(?:^|;[[:space:]]*)', key, '[[:space:]]+"([^"]*)"')
  z <- regmatches(x, regexec(pat, x, perl=TRUE))
  vapply(z, function(y) if (length(y) >= 2L) y[2L] else NA_character_, character(1))
}

# Search recursively because DESeq2 writes one folder per comparison.
result_files <- list.files(DTE_DIR, pattern='\\.tsv(\\.gz)?$', full.names=TRUE, recursive=TRUE)
result_files <- result_files[grepl(RESULT_PATTERN, basename(result_files))]
if (!length(result_files)) {
  message('No matching DTE files found; listing candidate TSV files:')
  print(head(list.files(DTE_DIR, pattern='\\.tsv(\\.gz)?$', recursive=TRUE), 40))
  stop('Adjust RESULT_PATTERN to match your actual DESeq2 result filenames.')
}
message('DTE result files found: ', length(result_files))

# Collect only IsoMatch transcript IDs needed for selected result files.
result_ids <- unique(unlist(lapply(result_files, function(f) {
  z <- fread(f, select='transcript_id', showProgress=FALSE)
  if (!'transcript_id' %in% names(z)) stop('Missing transcript_id in ', f)
  as.character(z$transcript_id)
}), use.names=FALSE))
message('Unique target transcripts: ', length(result_ids))

# Parse only needed IsoMatch transcript/exon features, streaming the very large GTF.
message('Reading target IsoMatch transcript/exon features...')
con <- gzfile(ISOMATCH_GTF, 'rt')
iso_parts <- list(); ni <- 0L
repeat {
  lines <- readLines(con, n=CHUNK_LINES, warn=FALSE)
  if (!length(lines)) break
  lines <- lines[!startsWith(lines, '#') & grepl('\\t(transcript|exon)\\t', lines)]
  if (!length(lines)) next
  f <- tstrsplit(lines, '\t', fixed=TRUE, keep=1:9)
  tid <- attr_get(f[[9]], 'transcript_id')
  sel <- !is.na(tid) & tid %chin% result_ids
  if (!any(sel)) next
  ni <- ni + 1L
  iso_parts[[ni]] <- data.table(transcript_id=tid[sel], chr=f[[1]][sel],
    feature=f[[3]][sel], start=as.integer(f[[4]][sel]), end=as.integer(f[[5]][sel]),
    strand=f[[7]][sel], isomatch_gene_id=attr_get(f[[9]][sel], 'gene_id'))
}
close(con)
if (!length(iso_parts)) stop('No requested transcript IDs found in IsoMatch GTF')
iso <- rbindlist(iso_parts); rm(iso_parts); gc()
# GTF files can contain malformed coordinates. Exclude them before constructing IRanges.
valid_iso <- !is.na(iso$start) & !is.na(iso$end) &
  iso$start >= 1L & iso$end >= iso$start &
  !is.na(iso$chr) & nzchar(iso$chr) & iso$strand %chin% c('+','-')
message('Excluded invalid IsoMatch transcript/exon records: ', sum(!valid_iso))
if (any(!valid_iso)) fwrite(iso[!valid_iso], file.path(OUT_DIR, 'invalid_IsoMatch_GTF_records.tsv'), sep='\t')
iso <- iso[valid_iso]
iso_tx <- unique(iso[feature=='transcript', .(transcript_id, isomatch_gene_id, chr, strand)])
iso_ex <- unique(iso[feature=='exon', .(transcript_id, chr, start, end, strand)])
if (!nrow(iso_ex)) stop('No IsoMatch exon records found')
setorder(iso_ex, transcript_id, start, end)

# Read GENCODE gene and exon records (not entire GTF into memory at once).
message('Reading GENCODE v38 genes and exons...')
con <- file(GENCODE_GTF, 'rt')
gene_parts <- list(); exon_parts <- list(); ng <- 0L; ne <- 0L
repeat {
  lines <- readLines(con, n=CHUNK_LINES, warn=FALSE)
  if (!length(lines)) break
  lines <- lines[!startsWith(lines, '#') & grepl('\\t(gene|exon)\\t', lines)]
  if (!length(lines)) next
  f <- tstrsplit(lines, '\t', fixed=TRUE, keep=1:9)
  gene_id <- attr_get(f[[9]], 'gene_id')
  feature <- f[[3]]
  is_gene <- feature=='gene' & !is.na(gene_id)
  if (any(is_gene)) {
    ng <- ng+1L
    gene_parts[[ng]] <- data.table(chr=f[[1]][is_gene], start=as.integer(f[[4]][is_gene]),
      end=as.integer(f[[5]][is_gene]), strand=f[[7]][is_gene],
      gencode_gene_id=gene_id[is_gene], gene_name=attr_get(f[[9]][is_gene], 'gene_name'))
  }
  is_exon <- feature=='exon' & !is.na(gene_id)
  if (any(is_exon)) {
    ne <- ne+1L
    exon_parts[[ne]] <- data.table(chr=f[[1]][is_exon], start=as.integer(f[[4]][is_exon]),
      end=as.integer(f[[5]][is_exon]), strand=f[[7]][is_exon],
      gencode_gene_id=gene_id[is_exon],
      gencode_transcript_id=attr_get(f[[9]][is_exon], 'transcript_id'))
  }
}
close(con)
genes <- unique(rbindlist(gene_parts), by='gencode_gene_id')
exons <- unique(rbindlist(exon_parts))
valid_gene <- !is.na(genes$start) & !is.na(genes$end) & genes$start >= 1L &
  genes$end >= genes$start & !is.na(genes$chr) & genes$strand %chin% c('+','-')
valid_exon <- !is.na(exons$start) & !is.na(exons$end) & exons$start >= 1L &
  exons$end >= exons$start & !is.na(exons$chr) & exons$strand %chin% c('+','-')
message('Excluded invalid GENCODE gene/exon records: ', sum(!valid_gene), '/', sum(!valid_exon))
genes <- genes[valid_gene]
exons <- exons[valid_exon & gencode_gene_id %chin% genes$gencode_gene_id]
rm(gene_parts, exon_parts); gc()
message('GENCODE genes: ', nrow(genes), '; exon records: ', nrow(exons))

# Precompute gene exon unions and annotated splice junction sets.
# Splitting by gene avoids counting shared bases repeatedly across transcript isoforms.
message('Building exon unions and annotated junctions...')
exons[, gene_index := match(gencode_gene_id, genes$gencode_gene_id)]
exons <- exons[!is.na(gene_index)]
exons[, transcript_key := paste(gencode_gene_id, gencode_transcript_id, sep='|')]
setorder(exons, transcript_key, start, end)
gj <- exons[, if (.N >= 2L) .(left=end[-.N], right=start[-1L]) else NULL,
            by=.(gene_index, transcript_key)]
gj <- unique(gj[, .(gene_index, left, right)])
gj[, jkey := paste(gene_index, left, right, sep=':')]
all_jkeys <- unique(gj$jkey)
rm(gj)

# Match only shared sequence levels to avoid Seqinfo mismatch warnings.
shared_chr <- intersect(unique(iso_ex$chr), unique(exons$chr))
message('Shared sequence levels: ', length(shared_chr))
iso_ex <- iso_ex[chr %chin% shared_chr]
# Candidate genes: same-strand overlap with the IsoMatch exonic span.
# Calculate exon-level intersection using unioned GENCODE exons for each gene.
ggr <- GRanges(seqnames=genes$chr, ranges=IRanges(genes$start, genes$end), strand=genes$strand)
iso_gr <- GRanges(seqnames=iso_ex$chr, ranges=IRanges(iso_ex$start, iso_ex$end), strand=iso_ex$strand)
# Avoid duplicating exon intervals by reducing GENCODE exons per gene.
exon_groups <- split(GRanges(seqnames=exons$chr, ranges=IRanges(exons$start, exons$end),
                             strand=exons$strand), exons$gene_index)
union_exons <- endoapply(exon_groups, reduce)
union_gr <- unlist(union_exons, use.names=FALSE)
union_gene <- rep(as.integer(names(union_exons)), lengths(union_exons))
rm(exon_groups); gc()

hits <- findOverlaps(iso_gr, union_gr, ignore.strand=FALSE)
q <- queryHits(hits); s <- subjectHits(hits)
inter <- pmax(0L, pmin(end(iso_gr)[q], end(union_gr)[s]) -
                   pmax(start(iso_gr)[q], start(union_gr)[s]) + 1L)
score <- data.table(transcript_id=iso_ex$transcript_id[q], gene_index=union_gene[s], overlap_bp=inter)
score <- score[, .(overlap_bp=sum(overlap_bp)), by=.(transcript_id, gene_index)]
iso_length <- iso_ex[, .(iso_exonic_bp=sum(end-start+1L)), by=transcript_id]
score <- merge(score, iso_length, by='transcript_id')
score[, exon_overlap_fraction := pmin(1, overlap_bp / iso_exonic_bp)]

# Calculate annotated intron-junction compatibility for each candidate gene.
iso_j <- iso_ex[, if (.N >= 2L) .(left=end[-.N], right=start[-1L]) else NULL, by=transcript_id]
iso_nj <- iso_j[, .(n_junctions=.N), by=transcript_id]
score <- merge(score, iso_nj, by='transcript_id', all.x=TRUE)
score[is.na(n_junctions), n_junctions := 0L]
if (nrow(iso_j)) {
  sj <- merge(score[, .(transcript_id, gene_index)], iso_j, by='transcript_id', allow.cartesian=TRUE)
  sj[, jkey := paste(gene_index, left, right, sep=':')]
  sj[, supported := jkey %chin% all_jkeys]
  ns <- sj[, .(n_supported_junctions=sum(supported)), by=.(transcript_id, gene_index)]
  score <- merge(score, ns, by=c('transcript_id','gene_index'), all.x=TRUE)
} else score[, n_supported_junctions := 0L]
score[is.na(n_supported_junctions), n_supported_junctions := 0L]
score[, junction_fraction := fifelse(n_junctions > 0, n_supported_junctions / n_junctions, NA_real_)]

# Conservative priority: all junctions supported, then exon overlap. If junctions
# are novel, allow a provisional unique exon-overlap assignment, explicitly labeled.
score[, tier := fifelse(n_junctions > 0 & n_supported_junctions == n_junctions &
                          exon_overlap_fraction >= MIN_EXON_OVERLAP, 2L,
                        fifelse(exon_overlap_fraction >= MIN_EXON_OVERLAP, 1L, 0L))]
setorder(score, transcript_id, -tier, -exon_overlap_fraction, gene_index)
ranked <- score[, {
  best <- .SD[1L]
  runner <- if (.N > 1L) .SD[2L] else NULL
  tie <- !is.null(runner) && runner$tier == best$tier &&
    (best$exon_overlap_fraction - runner$exon_overlap_fraction) < MIN_WINNER_MARGIN
  .(gene_index=best$gene_index, tier=best$tier,
    exon_overlap_fraction=best$exon_overlap_fraction,
    n_junctions=best$n_junctions,
    n_supported_junctions=best$n_supported_junctions,
    candidate_gene_count=.N, ambiguous=tie)
}, by=transcript_id]
ranked[, annotation_status := fifelse(tier==0L, 'insufficient_exon_support',
                              fifelse(ambiguous, 'ambiguous_gene_assignment',
                              fifelse(tier==2L, 'junction_supported_gene',
                                      'provisional_exon_supported_gene')))]
ranked[, `:=`(gencode_gene_id=genes$gencode_gene_id[gene_index],
               gene_name=genes$gene_name[gene_index])]
ranked[annotation_status %chin% c('insufficient_exon_support','ambiguous_gene_assignment'),
       `:=`(gencode_gene_id=NA_character_, gene_name=NA_character_)]

annotation <- merge(iso_tx, ranked[, .(transcript_id, gencode_gene_id, gene_name,
  annotation_status, exon_overlap_fraction, n_junctions, n_supported_junctions,
  candidate_gene_count)], by='transcript_id', all.x=TRUE)
annotation[is.na(annotation_status), annotation_status := 'no_same_strand_exon_overlap']
annotation <- unique(annotation, by='transcript_id')

# Distinguish strong junction evidence from provisional gene assignment.
# A full intron chain match to one GENCODE transcript is not established by
# individually supported junctions. Do not label these as exact transcript matches.
map_file <- file.path(OUT_DIR, 'IsoMatch_to_GENCODE_structure_aware.tsv.gz')
fwrite(annotation, map_file, sep='\t')
fwrite(annotation[, .N, by=annotation_status][order(-N)],
       file.path(OUT_DIR, 'annotation_status_summary.tsv'), sep='\t')
message('Saved annotation map: ', map_file)

for (f in result_files) {
  message('Annotating ', f)
  d <- fread(f)
  d[, .original_order := seq_len(.N)]
  d <- merge(d, annotation[, .(transcript_id, gencode_gene_id, gene_name,
    annotation_status, exon_overlap_fraction, n_junctions, n_supported_junctions,
    candidate_gene_count)], by='transcript_id', all.x=TRUE, sort=FALSE)
  # Preserve the original order of DESeq2 results.
  setorder(d, .original_order)
  d[, .original_order := NULL]
  relative <- substring(f, nchar(DTE_DIR)+2L)
  dest <- file.path(OUT_DIR, dirname(relative))
  dir.create(dest, recursive=TRUE, showWarnings=FALSE)
  outfile <- file.path(dest, paste0(sub('\\.tsv(\\.gz)?$', '', basename(f)), '_structure_annotated.tsv.gz'))
  fwrite(d, outfile, sep='\t')
}
message('Done. Output directory: ', OUT_DIR)
