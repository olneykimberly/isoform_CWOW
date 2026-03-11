# Read in annotation file
gtf.file <- paste0("/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/GRCh38/gencode.v38.annotation.sorted.gtf") 
gtf.gr <- rtracklayer::import(gtf.file)
# save gtf as data frame
gtf.df <- as.data.frame(gtf.gr)
# get gene id, transcript id, gene name, seqname which is chromosome, and biotype from gtf
genes <-
  gtf.df[, c("seqnames",
             "width",
             "gene_id",
             "gene_type",
             "gene_name",
             "transcript_id",
             "transcript_type", 
             "transcript_name")]
# output as a data frame
write.table(
  genes,
  "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/genes.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)