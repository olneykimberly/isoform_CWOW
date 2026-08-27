# ============================================================
# CWOW PacBio Kinnex cohort metadata summary
# n = 276
# ============================================================

library(tidyverse)
library(scales)
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts")
# ------------------------------------------------------------
# 1. Read metadata
# ------------------------------------------------------------

metadata <- read.csv(
  "CWOW_full_long_read_cohort_n276.csv",
  check.names = FALSE
)

# ------------------------------------------------------------
# 2. Clean disease and sex variables
# ------------------------------------------------------------

metadata <- metadata %>%
  mutate(
    Disease = factor(
      group,
      levels = c(
        "CONTROL",
        "AD",
        "LBD_S",
        "LBD_AS",
        "LBD_ATS"
      )
    ),

    # Gender currently contains Female/female and Male/male
    Sex = case_when(
      str_to_lower(Gender) == "female" ~ "Female",
      str_to_lower(Gender) == "male"   ~ "Male",
      TRUE ~ NA_character_
    ),

    Sex = factor(
      Sex,
      levels = c("Female", "Male")
    )
  )


# ------------------------------------------------------------
# 3. Sample counts by disease
# ------------------------------------------------------------

disease_counts <- metadata %>%
  count(Disease, name = "N")

print(disease_counts)


# ------------------------------------------------------------
# 4. Sample counts by disease and sex
# ------------------------------------------------------------

disease_sex_counts <- metadata %>%
  count(Disease, Sex, name = "N") %>%
  group_by(Disease) %>%
  mutate(
    Total = sum(N),
    Percent = N / Total * 100
  ) %>%
  ungroup()

print(disease_sex_counts)


# ------------------------------------------------------------
# 5. Publication-friendly summary table
# ------------------------------------------------------------

cohort_summary <- metadata %>%
  count(Disease, Sex) %>%
  pivot_wider(
    names_from = Sex,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    Total = Female + Male
  ) %>%
  select(
    Disease,
    Total,
    Female,
    Male
  )

print(cohort_summary)

write.csv(
  cohort_summary,
  "CWOW_Kinnex_cohort_summary.csv",
  row.names = FALSE
)


# ============================================================
# FIGURE 1
# Number of samples in each disease group
# ============================================================

p_disease <- ggplot(
  disease_counts,
  aes(x = Disease, y = N)
) +
  geom_col(
    width = 0.7,
    fill = "grey40"
  ) +
  geom_text(
    aes(label = N),
    vjust = -0.4,
    size = 5
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Number of individuals",
    title = "PacBio Kinnex long-read RNA-seq cohort"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

p_disease

ggsave(
  "CWOW_Kinnex_samples_by_disease.png",
  p_disease,
  width = 7,
  height = 5,
  dpi = 600
)


# ============================================================
# FIGURE 2
# Sex distribution within each disease group
# ============================================================

p_sex <- ggplot(
  disease_sex_counts,
  aes(
    x = Disease,
    y = N,
    fill = Sex
  )
) +
  geom_col(
    width = 0.7
  ) +
  geom_text(
    aes(label = N),
    position = position_stack(vjust = 0.5),
    size = 4.5
  ) +
  labs(
    x = NULL,
    y = "Number of individuals",
    fill = "Sex",
    title = "Sex distribution across neuropathological groups"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

p_sex

ggsave(
  "CWOW_Kinnex_sex_by_disease.png",
  p_sex,
  width = 7,
  height = 5,
  dpi = 600
)


# ============================================================
# FIGURE 3
# Optional: percentage composition
# ============================================================

p_sex_percent <- ggplot(
  disease_sex_counts,
  aes(
    x = Disease,
    y = N,
    fill = Sex
  )
) +
  geom_col(
    position = "fill",
    width = 0.7
  ) +
  geom_text(
    aes(
      label = paste0(round(Percent), "%")
    ),
    position = position_fill(vjust = 0.5),
    size = 4.5
  ) +
  scale_y_continuous(
    labels = percent_format()
  ) +
  labs(
    x = NULL,
    y = "Percent of disease group",
    fill = "Sex",
    title = "Sex composition across neuropathological groups"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

p_sex_percent

ggsave(
  "CWOW_Kinnex_sex_percent_by_disease.png",
  p_sex_percent,
  width = 7,
  height = 5,
  dpi = 600
)
