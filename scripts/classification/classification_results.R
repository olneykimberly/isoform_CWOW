# ============================================================
# CWOW Kinnex long-read transcriptome summary
# ============================================================

library(tidyverse)
library(readxl)
install.packages("janitor")
library(janitor)

# ------------------------------------------------------------
# USER INPUT
# ------------------------------------------------------------

# IMPORTANT:
# This should be the FINAL COHORT-LEVEL SQANTI3 classification
# table after your chosen high-confidence filtering.
#
# Do NOT simply concatenate all 276 sample classification
# files without collapsing transcript IDs across the cohort.

sqanti_file <- "CWOW_cohort_filtered_classification.txt"


# ------------------------------------------------------------
# 1. Read SQANTI3 classification table
# ------------------------------------------------------------

sqanti <- read.delim(
  sqanti_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# See available columns
names(sqanti)


# ------------------------------------------------------------
# 2. Identify relevant columns
# ------------------------------------------------------------

# Typical SQANTI3 columns include:
#
# isoform
# chrom
# strand
# length
# exons
# structural_category
# associated_gene
# associated_transcript
#
# Verify before proceeding.

sqanti_summary <- sqanti %>%
  filter(
    !is.na(isoform),
    !is.na(structural_category)
  )


# ------------------------------------------------------------
# 3. Total high-confidence transcripts
# ------------------------------------------------------------

N_transcripts <- n_distinct(
  sqanti_summary$isoform
)


# ------------------------------------------------------------
# 4. Total genes represented
# ------------------------------------------------------------

N_genes <- sqanti_summary %>%
  filter(
    !is.na(associated_gene),
    associated_gene != ""
  ) %>%
  distinct(associated_gene) %>%
  nrow()


# ------------------------------------------------------------
# 5. Structural categories
# ------------------------------------------------------------

category_summary <- sqanti_summary %>%
  count(
    structural_category,
    name = "N"
  ) %>%
  mutate(
    Percent = N / N_transcripts * 100
  ) %>%
  arrange(desc(N))

print(category_summary)


# ------------------------------------------------------------
# Full-splice match
# ------------------------------------------------------------

FSM <- category_summary %>%
  filter(
    structural_category == "full-splice_match"
  )

N_FSM <- FSM$N
PCT_FSM <- FSM$Percent


# ------------------------------------------------------------
# Novel-in-catalog
# ------------------------------------------------------------

NIC <- category_summary %>%
  filter(
    structural_category == "novel_in_catalog"
  )

N_NIC <- NIC$N
PCT_NIC <- NIC$Percent


# ------------------------------------------------------------
# Novel-not-in-catalog
# ------------------------------------------------------------

NNC <- category_summary %>%
  filter(
    structural_category == "novel_not_in_catalog"
  )

N_NNC <- NNC$N
PCT_NNC <- NNC$Percent


# ------------------------------------------------------------
# 6. Print abstract-ready sentence
# ------------------------------------------------------------

cat(
  paste0(
    "Across the cohort, we identified ",
    format(N_transcripts, big.mark = ","),
    " high-confidence transcripts from ",
    format(N_genes, big.mark = ","),
    " genes, including ",
    format(N_FSM, big.mark = ","),
    " (", round(PCT_FSM, 1), "%) full-splice matches, ",
    format(N_NIC, big.mark = ","),
    " (", round(PCT_NIC, 1), "%) novel-in-catalog transcripts, and ",
    format(N_NNC, big.mark = ","),
    " (", round(PCT_NNC, 1),
    "%) novel-not-in-catalog transcripts.\n"
  )
)