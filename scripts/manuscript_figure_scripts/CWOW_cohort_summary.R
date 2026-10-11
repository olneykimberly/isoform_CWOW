# ================================================================
# CWOW KINNEX LONG-READ COHORT: FINAL 274 DONORS
# Run interactively in Posit Workbench (Source).
# No command-line arguments required.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

# 1. PATHS AND SETTINGS -------------------------------------------
project_dir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW"
metadata_file <- file.path(project_dir, "scripts", "CWOW_full_long_read_cohort_n276.csv")
manifest_file <- file.path(project_dir, "scripts/isomatch_outputs/manifest_allSamples_deduplicated.txt")
output_dir <- file.path(project_dir, "results/cohort_demographics")
expected_n <- 274L

# Set metadata_file above to its actual location if different.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(metadata_file), file.exists(manifest_file))

group_levels <- c("CONTROL", "AD", "LBD_S", "LBD_AS", "LBD_ATS")
group_labels <- c("CONTROL" = "Control", "AD" = "AD(AT)",
                  "LBD_S" = "LBD(S)", "LBD_AS" = "LBD(AS)",
                  "LBD_ATS" = "LBD(ATS)")
# Consistent with Figure_QC_examples.R
group_colors <- c("CONTROL" = "#377EB8", "AD" = "#E41A1C",
                  "LBD_S" = "#6A3D9A", "LBD_AS" = "#9E9AC8",
                  "LBD_ATS" = "#636363")
sex_colors <- c("Female" = "#7B3294", "Male" = "#E69F00")

publication_theme <- function() {
  theme_classic(base_size = 13) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.position = "top", legend.title = element_text(face = "bold"))
}
save_figure <- function(plot, stem, width = 7.5, height = 5) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 600, bg = "white")
}

# 2. LOAD METADATA ------------------------------------------------
meta <- read_csv(metadata_file, col_types = cols(.default = col_character()),
                 show_col_types = FALSE)
required <- c("Study-Subject-ID", "NPID", "group", "Gender")
if (!all(required %in% names(meta))) {
  stop("Missing metadata columns: ", paste(setdiff(required, names(meta)), collapse = ", "))
}
meta <- meta %>% mutate(across(all_of(required), str_trim))
if (anyDuplicated(meta$NPID) || anyDuplicated(meta$`Study-Subject-ID`)) {
  stop("Duplicate metadata participant identifiers: inspect the source CSV.")
}

# Donor identifier: NA19-060 and NA19060-LBD-ATS both become NA19060.
# This is the sample identifier convention used in the cohort metadata.
extract_donor <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("NA([0-9]{2})-([0-9]{3})", "NA\\1\\2", x, perl = TRUE)
  str_extract(x, "NA[0-9]{5}")
}
meta <- meta %>% mutate(donor_id = extract_donor(NPID))
if (anyNA(meta$donor_id) || anyDuplicated(meta$donor_id)) {
  stop("Could not derive unique donor IDs from metadata NPID. Inspect NPID formatting.")
}

# 3. LOAD AND MATCH MANIFEST --------------------------------------
# Manifest may be a list of paths, IDs, or tab-separated records.
# We deliberately DO NOT use row order or randomly drop two samples.
manifest_lines <- readLines(manifest_file, warn = FALSE)
manifest_lines <- trimws(manifest_lines)
manifest_lines <- manifest_lines[nzchar(manifest_lines) &
                                   !grepl("^#", manifest_lines)]
if (!length(manifest_lines)) stop("Manifest is empty.")
cat("\nManifest preview:\n")
print(head(manifest_lines, 5))

# Extract ALL distinct donor IDs on each line to detect ambiguous lines.
extract_all_donors <- function(x) {
  x <- toupper(x)
  x <- gsub("NA([0-9]{2})-([0-9]{3})", "NA\\1\\2", x, perl = TRUE)
  unique(str_extract_all(x, "NA[0-9]{5}")[[1]])
}
line_donors <- lapply(manifest_lines, extract_all_donors)
manifest_audit <- tibble(
  manifest_line = seq_along(manifest_lines),
  record = manifest_lines,
  extracted_ids = vapply(line_donors, paste, collapse = ";", FUN.VALUE = character(1)),
  n_ids = lengths(line_donors)
)
write_csv(manifest_audit, file.path(output_dir, "CWOW_manifest_matching_audit.csv"))

if (any(manifest_audit$n_ids != 1L)) {
  print(manifest_audit %>% filter(n_ids != 1L), n = 20)
  stop("Manifest lines without exactly one donor ID. Review matching audit; ",
       "do not infer the cohort from an ambiguous manifest.")
}
manifest_ids <- unique(vapply(line_donors, `[[`, character(1), 1))
cat("Manifest lines:", length(manifest_lines),
    " | unique donors:", length(manifest_ids), "\n")
if (length(manifest_ids) != expected_n) {
  stop("Manifest identifies ", length(manifest_ids),
       " donors; expected ", expected_n, ". Check manifest content.")
}
unknown <- setdiff(manifest_ids, meta$donor_id)
if (length(unknown)) stop("Manifest donors missing in metadata: ",
                          paste(unknown, collapse = ", "))

final_meta <- meta %>% filter(donor_id %in% manifest_ids)
excluded_meta <- meta %>% filter(!donor_id %in% manifest_ids)
if (nrow(final_meta) != expected_n || n_distinct(final_meta$donor_id) != expected_n) {
  stop("Final cohort is not exactly 274 unique donors.")
}
cat("\nExcluded samples (metadata not in manifest):\n")
print(excluded_meta %>% select(`Study-Subject-ID`, NPID, group, Gender))
write_csv(excluded_meta, file.path(output_dir, "CWOW_excluded_samples.csv"))

# 4. CLEAN GROUP, SEX, AND NUMERIC VARIABLES ----------------------
final_meta <- final_meta %>% mutate(
  group = toupper(group),
  group = recode(group, "AD_AT" = "AD", "CTRL" = "CONTROL",
                 "LBD-S" = "LBD_S", "LBD-AS" = "LBD_AS",
                 "LBD-ATS" = "LBD_ATS"),
  Sex = case_when(tolower(Gender) %in% c("female", "f") ~ "Female",
                  tolower(Gender) %in% c("male", "m") ~ "Male",
                  TRUE ~ NA_character_)
)
if (anyNA(final_meta$Sex) || any(!final_meta$group %in% group_levels)) {
  stop("Unexpected or missing sex or disease group labels.")
}
final_meta <- final_meta %>% mutate(
  group = factor(group, levels = group_levels),
  Sex = factor(Sex, levels = c("Female", "Male"))
)

numeric_candidates <- c("Age", "PMI", "Braak.NFT", "Thal.amyloid", "Cing.LB")
for (nm in intersect(numeric_candidates, names(final_meta))) {
  final_meta[[nm]] <- suppressWarnings(as.numeric(final_meta[[nm]]))
}
write_csv(final_meta, file.path(output_dir, "CWOW_final_metadata_n274.csv"))

# 5. SAMPLE COUNTS AND TABLE 1 ------------------------------------
group_counts <- final_meta %>% count(group, name = "N", .drop = FALSE)
sex_counts <- final_meta %>% count(group, Sex, name = "N", .drop = FALSE) %>%
  group_by(group) %>% mutate(Group_N = sum(N), Percent = 100 * N / Group_N) %>%
  ungroup()

count_table <- sex_counts %>% select(group, Sex, N) %>%
  pivot_wider(names_from = Sex, values_from = N, values_fill = 0) %>%
  mutate(Total = Female + Male,
         `Female, n (%)` = sprintf("%d (%.1f%%)", Female, 100 * Female / Total),
         `Male, n (%)` = sprintf("%d (%.1f%%)", Male, 100 * Male / Total),
         Group = unname(group_labels[as.character(group)])) %>%
  select(Group, Total, `Female, n (%)`, `Male, n (%)`)
count_table <- bind_rows(count_table, tibble(
  Group = "Total", Total = nrow(final_meta),
  `Female, n (%)` = sprintf("%d (%.1f%%)", sum(final_meta$Sex == "Female"),
                            100 * mean(final_meta$Sex == "Female")),
  `Male, n (%)` = sprintf("%d (%.1f%%)", sum(final_meta$Sex == "Male"),
                          100 * mean(final_meta$Sex == "Male"))
))
print(count_table, width = Inf)
write_csv(count_table, file.path(output_dir, "Table1_sample_counts_by_sex.csv"))

fmt_iqr <- function(x, digits = 1) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  sprintf(paste0("%.", digits, "f [%.", digits, "f-%.", digits, "f]"),
          median(x), unname(quantile(x, .25)), unname(quantile(x, .75)))
}

# Only include available metadata variables; avoid assuming a column exists.
characteristic_vars <- intersect(numeric_candidates, names(final_meta))
make_characteristics <- function(data, label) {
  out <- data %>% summarise(
    Group = label, N = n(),
    `Female, n (%)` = sprintf("%d (%.1f%%)", sum(Sex == "Female"),
                              100 * mean(Sex == "Female")),
    `Male, n (%)` = sprintf("%d (%.1f%%)", sum(Sex == "Male"),
                            100 * mean(Sex == "Male")))
  for (nm in characteristic_vars) out[[paste0(nm, ", median [IQR]")]] <- fmt_iqr(data[[nm]])
  out
}
characteristics <- bind_rows(
  lapply(group_levels, function(g) {
    make_characteristics(final_meta %>% filter(group == g), group_labels[[g]])
  }),
  list(make_characteristics(final_meta, "Total"))
)
print(characteristics, width = Inf)
write_csv(characteristics, file.path(output_dir, "Table1_cohort_characteristics.csv"))

# 6. FIGURE 1: GROUP COUNTS ---------------------------------------
p_group <- ggplot(group_counts, aes(group, N, fill = group)) +
  geom_col(width = .72, color = "black", linewidth = .3) +
  geom_text(aes(label = N), vjust = -.45, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  labs(x = NULL, y = "Number of individuals",
       title = "PacBio Kinnex long-read RNA-seq cohort") +
  publication_theme() + theme(legend.position = "none")
print(p_group)
save_figure(p_group, "Figure1_samples_by_group")

# 7. FIGURE 2: SEX COUNTS -----------------------------------------
p_sex <- ggplot(sex_counts, aes(group, N, fill = Sex)) +
  geom_col(width = .72, color = "black", linewidth = .3) +
  geom_text(aes(label = ifelse(N > 0, N, "")),
            position = position_stack(vjust = .5),
            color = "white", fontface = "bold", size = 4.2) +
  geom_text(data = group_counts, aes(x = group, y = N, label = paste0("n=", N)),
            inherit.aes = FALSE, vjust = -.5, size = 3.7) +
  scale_fill_manual(values = sex_colors) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  labs(x = NULL, y = "Number of individuals", fill = "Sex",
       title = "Sex distribution across neuropathological groups") +
  publication_theme()
print(p_sex)
save_figure(p_sex, "Figure2_samples_by_group_and_sex")

# 8. FIGURE 3: PERCENTAGE BY SEX ----------------------------------
p_sex_percent <- ggplot(sex_counts, aes(group, Percent, fill = Sex)) +
  geom_col(width = .72, color = "black", linewidth = .3) +
  geom_text(aes(label = ifelse(N > 0, sprintf("%.1f%%", Percent), "")),
            position = position_stack(vjust = .5),
            color = "white", fontface = "bold", size = 3.6) +
  scale_fill_manual(values = sex_colors) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Individuals (%)", fill = "Sex",
       title = "Sex composition across neuropathological groups") +
  publication_theme()
print(p_sex_percent)
save_figure(p_sex_percent, "Figure3_sex_distribution_percent")

# 9. MANUSCRIPT DESCRIPTION AND FINAL VALIDATION ------------------
n_group <- setNames(group_counts$N, as.character(group_counts$group))
female_n <- sum(final_meta$Sex == "Female")
male_n <- sum(final_meta$Sex == "Male")
description <- sprintf(paste0(
  "The final long-read RNA sequencing cohort comprised %d individuals with ",
  "postmortem anterior cingulate cortex tissue. Neuropathological groups ",
  "included controls (n=%d), Alzheimer's disease with amyloid-beta and tau ",
  "pathology (AD(AT); n=%d), Lewy body disease with alpha-synuclein pathology ",
  "(LBD(S); n=%d), Lewy body disease with alpha-synuclein and amyloid-beta ",
  "pathology (LBD(AS); n=%d), and Lewy body disease with alpha-synuclein, ",
  "amyloid-beta, and tau pathology (LBD(ATS); n=%d). The cohort included ",
  "%d females (%.1f%%) and %d males (%.1f%%)."),
  nrow(final_meta), n_group[["CONTROL"]], n_group[["AD"]],
  n_group[["LBD_S"]], n_group[["LBD_AS"]], n_group[["LBD_ATS"]],
  female_n, 100 * female_n / nrow(final_meta),
  male_n, 100 * male_n / nrow(final_meta))

stopifnot(nrow(final_meta) == expected_n,
          n_distinct(final_meta$donor_id) == expected_n,
          sum(group_counts$N) == expected_n,
          sum(sex_counts$N) == expected_n,
          nrow(excluded_meta) == nrow(meta) - expected_n)
writeLines(c("CWOW LONG-READ COHORT SUMMARY", "", description, "",
             "Excluded metadata records:",
             capture.output(print(excluded_meta %>%
               select(`Study-Subject-ID`, NPID, group, Gender)))),
           file.path(output_dir, "CWOW_cohort_summary.txt"))
cat("\n", description, "\n\n", sep = "")
cat("Validated cohort: ", nrow(final_meta), " unique donors\n",
    "Female: ", female_n, " | Male: ", male_n, "\n",
    "Excluded: ", nrow(excluded_meta), "\n",
    "Results: ", output_dir, "\n", sep = "")

# 10. TABLE 1: PUBLICATION FORMAT (MEAN +/- SD AND PATHOLOGY BINS) ----
# Matches the columns in the example manuscript table.
# Counts reflect ONLY the manifest-verified 274 individuals.
# '-' indicates zero samples in a category; missing values are not
# assigned to a bin. The completeness audit is saved separately.

fmt_mean_sd <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  if (length(x) == 1L) return(sprintf(paste0("%.", digits, "f +/- NA"), x))
  sprintf(paste0("%.", digits, "f +/- %.", digits, "f"), mean(x), sd(x))
}
fmt_bin <- function(x) ifelse(x == 0L, "-", as.character(x))

make_paper_row <- function(d, group_name) {
  thal <- d$Thal.amyloid
  braak <- d$Braak.NFT
  tibble(
    Type = group_name,
    `No.` = nrow(d),
    `% male` = sprintf("%.1f%%", 100 * mean(d$Sex == "Male")),
    `Age, years mean +/- SD` = fmt_mean_sd(d$Age),
    `Cingulate LB (S) mean +/- SD` = fmt_mean_sd(d$Cing.LB),
    `Thal 0 to 1` = sum(thal %in% 0:1, na.rm = TRUE),
    `Thal 2 to 3` = sum(thal %in% 2:3, na.rm = TRUE),
    `Thal 4 to 5` = sum(thal %in% 4:5, na.rm = TRUE),
    `Braak 0 to III` = sum(braak %in% 0:3, na.rm = TRUE),
    `Braak IV to V` = sum(braak %in% 4:5, na.rm = TRUE),
    `Braak VI` = sum(braak == 6, na.rm = TRUE)
  )
}

paper_table_numeric <- bind_rows(lapply(group_levels, function(g) {
  make_paper_row(final_meta %>% filter(group == g), unname(group_labels[[g]]))
}))

# Confirm no pathology values fall outside the displayed bins.
pathology_audit <- final_meta %>% group_by(group) %>% summarise(
  N = n(),
  Thal_missing = sum(is.na(Thal.amyloid)),
  Thal_outside_0_5 = sum(!is.na(Thal.amyloid) &
                            !Thal.amyloid %in% 0:5),
  Braak_missing = sum(is.na(Braak.NFT)),
  Braak_outside_0_6 = sum(!is.na(Braak.NFT) &
                             !Braak.NFT %in% 0:6),
  Age_missing = sum(is.na(Age)),
  Cing_LB_missing = sum(is.na(Cing.LB)),
  .groups = "drop"
)
write_csv(pathology_audit, file.path(output_dir, "Table1_pathology_completeness_audit.csv"))
if (any(pathology_audit$Thal_outside_0_5 > 0 |
        pathology_audit$Braak_outside_0_6 > 0)) {
  warning("Some pathology values fall outside Table 1 bins; inspect completeness audit.")
}

paper_table <- paper_table_numeric %>% mutate(
  across(c(`Thal 0 to 1`, `Thal 2 to 3`, `Thal 4 to 5`,
           `Braak 0 to III`, `Braak IV to V`, `Braak VI`), fmt_bin)
)
print(paper_table, width = Inf)
write_csv(paper_table, file.path(output_dir, "Table1_publication_style_n274.csv"))

# Optional polished HTML table, if gt is installed.
if (requireNamespace("gt", quietly = TRUE)) {
  table_gt <- gt::gt(paper_table) %>%
    gt::tab_header(
      title = gt::md("**Table 1.** Anterior cingulate cortex sample demographics and pathological features"),
      subtitle = paste0("Total n = ", nrow(final_meta))
    ) %>%
    gt::tab_spanner(label = "Thal amyloid phase (A)",
                    columns = c(`Thal 0 to 1`, `Thal 2 to 3`, `Thal 4 to 5`)) %>%
    gt::tab_spanner(label = "Braak stage (T)",
                    columns = c(`Braak 0 to III`, `Braak IV to V`, `Braak VI`)) %>%
    gt::cols_label(
      `Age, years mean +/- SD` = "Age, years (mean +/- SD)",
      `Cingulate LB (S) mean +/- SD` = "Cingulate LB (mean +/- SD)",
      `Thal 0 to 1` = "0 to 1", `Thal 2 to 3` = "2 to 3",
      `Thal 4 to 5` = "4 to 5", `Braak 0 to III` = "0 to III",
      `Braak IV to V` = "IV to V", `Braak VI` = "VI"
    )
  gt::gtsave(table_gt, file.path(output_dir, "Table1_publication_style_n274.html"))
}
cat("\nPublication-style Table 1 saved. Review completeness audit for missing pathology values.\n")

# 10. RIN DISTRIBUTION BY NEUROPATHOLOGICAL GROUP ------------------
# Input: original RNA sample spreadsheet; only final 274 donors are plotted.
# Requires readxl (install once if needed: install.packages("readxl")).

rin_file <- file.path(project_dir, "scripts", "CWOW_Original_RNA_Samples.xlsx")
# Change rin_file above if the Excel workbook is stored elsewhere.
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required for the RIN plot. Install it with install.packages('readxl').")
}
if (!file.exists(rin_file)) {
  stop("RIN workbook not found: ", rin_file, "\nUpdate rin_file in Section 10.")
}

rin_raw <- readxl::read_excel(rin_file, sheet = 1)
rin_required <- c("DCL Sample ID", "DCL Patient ID", "Study Sample ID", "RIN")
if (!all(rin_required %in% names(rin_raw))) {
  stop("Missing required RIN columns: ",
       paste(setdiff(rin_required, names(rin_raw)), collapse = ", "))
}

# The study ID can be NA20-121 or NA20121; the patient ID can be
# CWOW_NA20121 or CWOW_NA20121-LBD-S. All normalize to NA20121.
rin_raw <- rin_raw %>%
  mutate(
    rin_donor_study = extract_donor(`Study Sample ID`),
    rin_donor_patient = extract_donor(`DCL Patient ID`),
    RIN = suppressWarnings(as.numeric(`RIN`)),
    donor_id = coalesce(rin_donor_study, rin_donor_patient)
  )

# Detect disagreements rather than arbitrarily selecting an identifier.
rin_disagree <- rin_raw %>%
  filter(!is.na(rin_donor_study), !is.na(rin_donor_patient),
         rin_donor_study != rin_donor_patient)
if (nrow(rin_disagree) > 0) {
  readr::write_csv(rin_disagree,
                   file.path(output_dir, "RIN_conflicting_sample_identifiers.csv"))
  stop("RIN workbook has conflicting study/patient IDs; inspect RIN_conflicting_sample_identifiers.csv")
}

rin_final_rows <- rin_raw %>% filter(donor_id %in% final_meta$donor_id)
# Multiple rows for the same donor must be reconciled manually.
rin_duplicates <- rin_final_rows %>% count(donor_id) %>% filter(n > 1)
if (nrow(rin_duplicates) > 0) {
  readr::write_csv(rin_final_rows %>% semi_join(rin_duplicates, by = "donor_id"),
                   file.path(output_dir, "RIN_duplicate_donor_rows.csv"))
  stop("Multiple RIN records for final-cohort donors; inspect RIN_duplicate_donor_rows.csv")
}

rin_data <- final_meta %>%
  select(donor_id, `Study-Subject-ID`, NPID, group, Sex) %>%
  left_join(
    rin_final_rows %>%
      select(donor_id, `DCL Sample ID`, `DCL Patient ID`, `Study Sample ID`, RIN),
    by = "donor_id"
  )

stopifnot(nrow(rin_data) == expected_n,
          n_distinct(rin_data$donor_id) == expected_n)

rin_missing <- rin_data %>% filter(is.na(RIN))
readr::write_csv(rin_data, file.path(output_dir, "RIN_final_cohort_n274.csv"))
readr::write_csv(rin_missing, file.path(output_dir, "RIN_missing_or_unmatched_donors.csv"))

rin_summary <- rin_data %>%
  group_by(group) %>%
  summarise(
    N_cohort = n(),
    N_with_RIN = sum(!is.na(RIN)),
    N_missing_RIN = sum(is.na(RIN)),
    Mean_RIN = if (all(is.na(RIN))) NA_real_ else mean(RIN, na.rm = TRUE),
    SD_RIN = if (sum(!is.na(RIN)) < 2) NA_real_ else sd(RIN, na.rm = TRUE),
    Median_RIN = if (all(is.na(RIN))) NA_real_ else median(RIN, na.rm = TRUE),
    Q1_RIN = if (all(is.na(RIN))) NA_real_ else as.numeric(quantile(RIN, 0.25, na.rm = TRUE)),
    Q3_RIN = if (all(is.na(RIN))) NA_real_ else as.numeric(quantile(RIN, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

print(rin_summary, width = Inf)
cat("\nRIN values available for", sum(!is.na(rin_data$RIN)),
    "of", nrow(rin_data), "final-cohort donors.\n")
if (nrow(rin_missing) > 0) {
  warning(nrow(rin_missing),
          " final-cohort donors lack matched RIN values. Review RIN_missing_or_unmatched_donors.csv")
}
readr::write_csv(rin_summary, file.path(output_dir, "Table_RIN_by_group.csv"))

rin_plot_data <- rin_data %>% filter(!is.na(RIN))
if (nrow(rin_plot_data) == 0) stop("No final-cohort RIN values available to plot.")

p_rin <- ggplot(rin_plot_data, aes(x = group, y = RIN, fill = group)) +
  geom_boxplot(width = 0.60, outlier.shape = NA,
               color = "black", linewidth = 0.55, alpha = 0.85) +
  geom_jitter(aes(color = group), width = 0.14,
              height = 0, size = 1.7, alpha = 0.65) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  labs(x = NULL, y = "RNA integrity number (RIN)",
       title = "RNA integrity by neuropathological group") +
  publication_theme() +
  theme(legend.position = "none")

print(p_rin)
save_figure(p_rin, "Figure_RIN_by_group", width = 7.5, height = 5)
cat("RIN figure, summary, and matching audits saved to: ", output_dir, "\n")

# 11. FINAL-COHORT RIN SUMMARY (OVERALL + BY GROUP) --------------
# Uses rin_data created in Section 10, not a nonexistent
# 'final_metadata' object. Only the 274 manifest-retained donors
# are included; donors without RIN are counted but excluded from
# RIN mean/median/SD calculations.
stopifnot(nrow(rin_data) == 274L,
          dplyr::n_distinct(rin_data$donor_id) == 274L)

rin_stats <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  valid <- x[is.finite(x)]
  tibble::tibble(
    N_with_RIN = length(valid),
    N_missing_RIN = length(x) - length(valid),
    Mean_RIN = if (length(valid)) mean(valid) else NA_real_,
    SD_RIN = if (length(valid) >= 2) stats::sd(valid) else NA_real_,
    Median_RIN = if (length(valid)) stats::median(valid) else NA_real_,
    Min_RIN = if (length(valid)) min(valid) else NA_real_,
    Max_RIN = if (length(valid)) max(valid) else NA_real_
  )
}

rin_by_group_final <- rin_data %>%
  group_by(group) %>%
  summarise(N_cohort = n(), rin_stats(RIN), .groups = "drop") %>%
  mutate(group = as.character(group))

rin_overall_final <- rin_data %>%
  summarise(N_cohort = n(), rin_stats(RIN)) %>%
  mutate(group = "Overall", .before = 1)

rin_report_final <- bind_rows(rin_by_group_final, rin_overall_final) %>%
  mutate(across(c(Mean_RIN, SD_RIN, Median_RIN, Min_RIN, Max_RIN),
                ~ round(.x, 2)))

cat("\nFINAL 274-SAMPLE COHORT: RIN REPORT\n")
print(rin_report_final, width = Inf)
readr::write_csv(rin_report_final,
                 file.path(output_dir, "Table_RIN_final_n274.csv"))
cat("\nMissing RIN donor details: ",
    file.path(output_dir, "RIN_missing_or_unmatched_donors.csv"), "\n", sep = "")
