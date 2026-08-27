setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts")
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

outdir <- "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/isomatch_cohort"

# qc should already exist from the previous step:
# classification + support joined by transcript_id

# ------------------------------------------------------------
# 1. Verify basic QC
# ------------------------------------------------------------

cat("Total transcripts:", format(nrow(qc), big.mark = ","), "\n")
cat("Missing sample support:", sum(is.na(qc$n_samples)), "\n")

table(qc$structural_category)


# ------------------------------------------------------------
# 2. Define recurrence thresholds
# ------------------------------------------------------------

thresholds <- c(1, 2, 3, 5, 10, 20, 50)


# ------------------------------------------------------------
# 3. Count retained transcripts at each threshold
#    by structural category
# ------------------------------------------------------------

filter_sensitivity <- lapply(
    thresholds,
    function(k) {

        qc %>%
            group_by(structural_category) %>%
            summarise(
                retained = sum(n_samples >= k),
                total = n(),
                pct_retained = 100 * retained / total,
                .groups = "drop"
            ) %>%
            mutate(
                min_samples = k
            )
    }
) %>%
    bind_rows()


filter_sensitivity <- filter_sensitivity %>%
    select(
        min_samples,
        structural_category,
        retained,
        total,
        pct_retained
    ) %>%
    arrange(
        min_samples,
        structural_category
    )

print(filter_sensitivity)

fwrite(
    filter_sensitivity,
    file.path(
        outdir,
        "CWOW_filter_sensitivity_by_category.tsv"
    ),
    sep = "\t"
)

main_categories <- c(
    "full-splice_match",
    "incomplete-splice_match",
    "novel_in_catalog",
    "novel_not_in_catalog"
)

main_filter_sensitivity <- filter_sensitivity %>%
    filter(
        structural_category %in% main_categories
    )

print(main_filter_sensitivity)

overall_sensitivity <- tibble(
    min_samples = thresholds
) %>%
    rowwise() %>%
    mutate(
        N_transcripts = sum(qc$n_samples >= min_samples),
        pct_of_catalog = 100 * N_transcripts / nrow(qc)
    ) %>%
    ungroup()

print(overall_sensitivity)

fwrite(
    overall_sensitivity,
    file.path(
        outdir,
        "CWOW_filter_sensitivity_overall.tsv"
    ),
    sep = "\t"
)

p_support <- main_filter_sensitivity %>%
    ggplot(
        aes(
            x = min_samples,
            y = retained,
            group = structural_category
        )
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    scale_x_log10(
        breaks = thresholds
    ) +
    scale_y_continuous(
        labels = comma
    ) +
    facet_wrap(
        ~ structural_category,
        scales = "free_y",
        ncol = 2
    ) +
    labs(
        x = "Minimum number of supporting individuals",
        y = "Retained transcripts",
        title = "Effect of cohort recurrence filtering",
        subtitle = "CWOW IsoMatch merged transcriptome"
    ) +
    theme_classic(base_size = 13)

p_support

ggsave(
    file.path(
        outdir,
        "CWOW_filter_sensitivity.png"
    ),
    p_support,
    width = 9,
    height = 7,
    dpi = 300
)

singleton_summary <- qc %>%
    group_by(structural_category) %>%
    summarise(
        total = n(),
        singleton = sum(n_samples == 1),
        pct_singleton = 100 * singleton / total,
        median_support = median(n_samples),
        mean_support = mean(n_samples),
        .groups = "drop"
    ) %>%
    arrange(desc(pct_singleton))

print(singleton_summary)

print(
    main_filter_sensitivity,
    n = Inf
)

# ============================================================
# Compare candidate definitions of the cohort transcriptome
# ============================================================

candidate_summary <- bind_rows(

    # --------------------------------------------------------
    # Candidate A:
    # FSM retained regardless of recurrence
    # NIC/NNC require >=2 individuals
    # --------------------------------------------------------

    qc %>%
        filter(
            structural_category == "full-splice_match" |
            (
                structural_category %in%
                    c(
                        "novel_in_catalog",
                        "novel_not_in_catalog"
                    ) &
                n_samples >= 2
            )
        ) %>%
        mutate(candidate = "FSM any; novel >=2"),

    # --------------------------------------------------------
    # Candidate B:
    # FSM retained regardless of recurrence
    # NIC/NNC require >=3 individuals
    # --------------------------------------------------------

    qc %>%
        filter(
            structural_category == "full-splice_match" |
            (
                structural_category %in%
                    c(
                        "novel_in_catalog",
                        "novel_not_in_catalog"
                    ) &
                n_samples >= 3
            )
        ) %>%
        mutate(candidate = "FSM any; novel >=3"),

    # --------------------------------------------------------
    # Candidate C:
    # FSM retained regardless of recurrence
    # NIC/NNC require >=5 individuals
    # --------------------------------------------------------

    qc %>%
        filter(
            structural_category == "full-splice_match" |
            (
                structural_category %in%
                    c(
                        "novel_in_catalog",
                        "novel_not_in_catalog"
                    ) &
                n_samples >= 5
            )
        ) %>%
        mutate(candidate = "FSM any; novel >=5")
) %>%

    group_by(
        candidate,
        structural_category
    ) %>%

    summarise(
        N = n(),
        .groups = "drop"
    ) %>%

    group_by(candidate) %>%

    mutate(
        total_catalog = sum(N),
        percent = 100 * N / total_catalog
    ) %>%

    ungroup()

print(candidate_summary, n = Inf)

hc_summary <- hc_ids[
  ,
  .(
    N = .N
  ),
  by = structural_category
][
  ,
  percent := 100 * N / sum(N)
][
  order(-N)
]

print(hc_summary)