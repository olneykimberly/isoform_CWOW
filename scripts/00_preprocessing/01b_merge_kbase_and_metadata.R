kbase <- read.delim("metadata/simple_kbase_name.tmp", header = FALSE)
meta <- read.delim("metadata/NPID_and_sex_chr.tmp")

merged_data <- merge(
  x = kbase,
  y = meta,
  by.x = "V2", 
  by.y = "NPID", 
  all = FALSE 
)

write.table(merged_data, "metadata/kbase_and_metadata.tsv", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
