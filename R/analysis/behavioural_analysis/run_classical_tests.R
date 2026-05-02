# =============================================================================
# run_classical_tests.R
# -----------------------------------------------------------------------------
# Classical inferential tests complementing the GAMM pipeline.
# Each test answers a question the GAMMs cannot address directly.
#
# Test 1  Paired t-test     ERD significance (pre vs post alpha power)
# Test 2  Paired t-test     Channel vs source GAMM fit quality (per-subject R²)
# Test 3  rmANOVA           BI_pre / CoG_pre / delta_ERD across ROIs
# Test 4  One-way ANOVA     Source GAMM deviance explained across ROIs
# Test 5  Rayleigh test     Phase uniformity for slow_phase and slow_phase_post
#
# Outputs  (<out_dir>/)
# ──────────────────────
#   test1_erd_ttest.csv              t / df / p / Cohen's d per ROI × band
#   test1_erd_plot.png               paired boxplot pre vs post per ROI
#   test2_fit_ttest.csv              t / df / p / Cohen's d per model concept
#   test2_fit_plot.png               paired dot-and-line plot
#   test3_rmanova_<metric>.csv       ANOVA table + pairwise post-hoc per metric
#   test3_rmanova_plot.png           faceted boxplot ROI × metric
#   test4_roi_anova_<model>.csv      ANOVA table + post-hoc per model concept
#   test4_roi_anova_plot.png         bar chart of mean R² per ROI
#   test5_rayleigh_<phase_col>.csv   z / p / R_mean per subject × ROI
#   test5_rayleigh_summary.csv       ROI-level summary (n_sig subjects)
#   test5_rayleigh_plot_<col>.png    mean resultant length heatmap + sig overlay
#   classical_tests_report.txt       plain-text narrative summary
# =============================================================================

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(purrr)
    library(tibble)
    library(ggplot2)
    library(ggpubr)    # for stat_compare_means and significance brackets
    library(broom)     # for tidy() on aov objects (Test 4)
})

# afex is preferred for rmANOVA (Greenhouse-Geisser correction).
# Falls back to base aov() with Error() if not installed.
.afex_ok <- requireNamespace("afex", quietly = TRUE)
if (.afex_ok) library(afex)

# =============================================================================
# USER SETTINGS
# =============================================================================
source_master_file  <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/source_pain_master.csv"
channel_master_file <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/alpha_pain_master.csv"
channel_fitted_file <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs_v2/trial_level_fitted_values_v2.csv"
source_out_dir      <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs_source"
out_dir             <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/classical_tests"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Which source metrics to test in rmANOVA (Test 3)
rmanova_metrics <- c("BI_pre", "CoG_pre", "delta_ERD")

# Which source model fitted CSVs to use for per-subject R² (Test 4)
# These match the model names produced by run_gamm_source.R
source_models_for_r2 <- c("m01_bi_pre", "m02_lr_pre", "m03_cog_pre",
                          "m07_delta_erd", "m08_n2p2_amp")

# Matched concept pairs for Test 2 (channel model name → source model name)
# Must match names in run_gamm_alpha_metrics_v2.R and run_gamm_source.R
concept_pairs <- tribble(
    ~concept,          ~channel_model,            ~source_model,
    "Log ratio",       "m02_sf_logratio",          "m02_lr_pre",
    "Balance index",   "m03_sf_balance",            "m01_bi_pre",
    "Slow fraction",   "m04_slow_frac",             "m03_cog_pre",
    "Tensor slow×fast","m05_tensor_slow_fast",       "m05_tensor_slow_fast"
)

alpha_level <- 0.05   # significance threshold for all tests

# =============================================================================
# HELPERS
# =============================================================================
clean_names_local <- function(x) {
    x %>% str_trim() %>%
        str_replace_all("\\s+", "_") %>%
        str_replace_all("\\^2",  "r2") %>%
        str_replace_all("\\.",   "_")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Cohen's d for paired differences
cohens_d_paired <- function(x, y) {
    diff  <- x - y
    mean(diff, na.rm = TRUE) / sd(diff, na.rm = TRUE)
}

# Compact paired t-test wrapper: returns one-row tibble
run_paired_ttest <- function(x, y, label) {
    complete <- !is.na(x) & !is.na(y)
    x <- x[complete]; y <- y[complete]
    if (length(x) < 3L) {
        return(tibble(label = label, n = length(x), t = NA, df = NA,
                      p_value = NA, cohens_d = NA, mean_diff = NA, se_diff = NA))
    }
    tt <- t.test(x, y, paired = TRUE)
    tibble(
        label     = label,
        n         = length(x),
        t         = round(tt$statistic, 3),
        df        = round(tt$parameter, 1),
        p_value   = tt$p.value,
        cohens_d  = round(cohens_d_paired(x, y), 3),
        mean_diff = round(tt$estimate, 4),
        se_diff   = round(diff(tt$conf.int) / (2 * 1.96), 4)
    )
}

# Per-subject R²  from a trial-level data frame with observed and fitted columns
per_subject_r2 <- function(df, obs_col, fit_col, id_col = "subjid_uid") {
    df %>%
        group_by(across(all_of(id_col))) %>%
        summarise(
            r2 = {
                o <- .data[[obs_col]]; f <- .data[[fit_col]]
                ok <- !is.na(o) & !is.na(f)
                if (sum(ok) < 3L) NA_real_ else cor(o[ok], f[ok])^2
            },
            .groups = "drop"
        )
}

# FDR-corrected pairwise t-tests across groups (long format)
pairwise_fdr <- function(df, value_col, group_col, id_col) {
    groups <- sort(unique(df[[group_col]]))
    pairs  <- combn(groups, 2, simplify = FALSE)
    
    map_dfr(pairs, function(p) {
        d1 <- df %>% filter(.data[[group_col]] == p[1])
        d2 <- df %>% filter(.data[[group_col]] == p[2])
        joined <- inner_join(
            d1 %>% select(all_of(c(id_col, value_col))),
            d2 %>% select(all_of(c(id_col, value_col))),
            by = id_col, suffix = c("_1", "_2")
        )
        x <- joined[[paste0(value_col, "_1")]]
        y <- joined[[paste0(value_col, "_2")]]
        run_paired_ttest(x, y, paste(p[1], "vs", p[2])) %>%
            mutate(group1 = p[1], group2 = p[2])
    }) %>%
        mutate(p_adj = p.adjust(p_value, method = "fdr"),
               sig   = p_adj < alpha_level)
}

# Significance label from p-value
sig_label <- function(p) {
    case_when(
        is.na(p)    ~ "n/a",
        p < 0.001   ~ "***",
        p < 0.01    ~ "**",
        p < 0.05    ~ "*",
        TRUE        ~ "ns"
    )
}

# =============================================================================
# LOAD DATA
# =============================================================================
message("Loading data...")

read_safe <- function(f) {
    if (!file.exists(f)) { message("[WARN] Not found: ", f); return(NULL) }
    df <- read_csv(f, show_col_types = FALSE)
    names(df) <- clean_names_local(names(df))
    df
}

src_master  <- read_safe(source_master_file)
chan_master  <- read_safe(channel_master_file)
chan_fitted  <- read_safe(channel_fitted_file)

if (is.null(src_master)) stop("source_pain_master.csv not found.")

all_rois <- sort(unique(src_master$roi))
message("ROIs: ", paste(all_rois, collapse = ", "))

# Subject-level GA per ROI (used by Tests 1, 3, 4)
src_subject_ga <- src_master %>%
    group_by(experiment_id, subjid_uid, roi) %>%
    summarise(
        across(
            any_of(c("pow_slow", "pow_fast", "pow_alpha",
                     "pow_slow_post", "pow_fast_post", "pow_alpha_post",
                     "BI_pre", "LR_pre", "CoG_pre", "psi_cog",
                     "delta_ERD", "ERD_slow", "ERD_fast",
                     "n2p2_amp", "n2_amp", "p2_amp",
                     "slow_phase", "slow_phase_post")),
            ~mean(.x, na.rm = TRUE)
        ),
        n_trials = n(),
        .groups  = "drop"
    )

# =============================================================================
# TEST 1: PAIRED t-TEST  —  ERD SIGNIFICANCE (pre vs post alpha power)
# =============================================================================
message("\n", strrep("=", 60))
message("TEST 1: ERD significance (pre vs post alpha power)")
message(strrep("=", 60))

test1_rows <- list()

for (this_roi in all_rois) {
    roi_df <- src_subject_ga %>% filter(roi == this_roi)
    
    # Slow alpha: pre vs post
    if (all(c("pow_slow", "pow_slow_post") %in% names(roi_df))) {
        test1_rows[[length(test1_rows) + 1L]] <-
            run_paired_ttest(roi_df$pow_slow, roi_df$pow_slow_post,
                             paste0(this_roi, " | slow_alpha")) %>%
            mutate(roi = this_roi, band = "slow")
    }
    # Fast alpha: pre vs post
    if (all(c("pow_fast", "pow_fast_post") %in% names(roi_df))) {
        test1_rows[[length(test1_rows) + 1L]] <-
            run_paired_ttest(roi_df$pow_fast, roi_df$pow_fast_post,
                             paste0(this_roi, " | fast_alpha")) %>%
            mutate(roi = this_roi, band = "fast")
    }
    # Total alpha: pre vs post
    if (all(c("pow_alpha", "pow_alpha_post") %in% names(roi_df))) {
        test1_rows[[length(test1_rows) + 1L]] <-
            run_paired_ttest(roi_df$pow_alpha, roi_df$pow_alpha_post,
                             paste0(this_roi, " | total_alpha")) %>%
            mutate(roi = this_roi, band = "total")
    }
}

test1_results <- bind_rows(test1_rows) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr"),
           sig   = sig_label(p_adj)) %>%
    arrange(roi, band)

write_csv(test1_results, file.path(out_dir, "test1_erd_ttest.csv"))
message("  Saved test1_erd_ttest.csv  (", nrow(test1_results), " rows)")

# ── Test 1 plot: paired boxplot pre vs post ──────────────────────────────────
test1_long <- src_master %>%
    select(subjid_uid, roi, experiment_id,
           any_of(c("pow_slow", "pow_slow_post",
                    "pow_fast", "pow_fast_post"))) %>%
    group_by(subjid_uid, roi, experiment_id) %>%
    summarise(across(everything(), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(
        cols      = any_of(c("pow_slow", "pow_slow_post",
                             "pow_fast", "pow_fast_post")),
        names_to  = "measure",
        values_to = "power"
    ) %>%
    mutate(
        band   = if_else(str_detect(measure, "slow"), "Slow alpha", "Fast alpha"),
        window = if_else(str_detect(measure, "post"), "Post-stim", "Pre-stim"),
        window = factor(window, levels = c("Pre-stim", "Post-stim"))
    ) %>%
    filter(!is.na(power))

if (nrow(test1_long) > 0) {
    p1 <- ggplot(test1_long, aes(x = window, y = power, fill = window)) +
        geom_boxplot(outlier.size = 0.5, width = 0.55, alpha = 0.8) +
        geom_line(aes(group = subjid_uid), alpha = 0.2, linewidth = 0.3, colour = "grey40") +
        stat_compare_means(method = "t.test", paired = TRUE,
                           label = "p.signif", size = 3.5,
                           comparisons = list(c("Pre-stim", "Post-stim"))) +
        facet_grid(band ~ roi, scales = "free_y") +
        scale_fill_manual(values = c("Pre-stim" = "steelblue", "Post-stim" = "tomato")) +
        labs(title = "Test 1: Pre vs Post Stimulus Alpha Power by ROI",
             subtitle = "Paired t-test; * p<0.05, ** p<0.01, *** p<0.001 (FDR corrected)",
             x = NULL, y = "Alpha power (source units)", fill = NULL) +
        theme_minimal(base_size = 10) +
        theme(legend.position = "top",
              axis.text.x = element_text(angle = 25, hjust = 1))
    
    ggsave(file.path(out_dir, "test1_erd_plot.png"),
           p1, width = max(8, length(all_rois) * 2), height = 6, dpi = 200)
    message("  Saved test1_erd_plot.png")
}

# =============================================================================
# TEST 2: PAIRED t-TEST  —  CHANNEL vs SOURCE FIT QUALITY (per-subject R²)
# =============================================================================
message("\n", strrep("=", 60))
message("TEST 2: Channel vs source GAMM fit quality (per-subject R²)")
message(strrep("=", 60))

test2_rows <- list()
test2_plot_data <- list()

for (i in seq_len(nrow(concept_pairs))) {
    concept    <- concept_pairs$concept[i]
    ch_model   <- concept_pairs$channel_model[i]
    src_model  <- concept_pairs$source_model[i]
    
    # ── Channel per-subject R² ─────────────────────────────────────────────────
    ch_r2 <- NULL
    fit_col_ch <- paste0(ch_model, "_fitted")
    if (!is.null(chan_fitted) && fit_col_ch %in% names(chan_fitted) &&
        "pain_rating" %in% names(chan_fitted) && "subjid_uid" %in% names(chan_fitted)) {
        ch_r2 <- per_subject_r2(chan_fitted, "pain_rating", fit_col_ch) %>%
            rename(r2_channel = r2)
    }
    
    if (is.null(ch_r2)) {
        message("  [SKIP] Channel fitted values not available for: ", ch_model)
        next
    }
    
    # ── Source per-subject R²: best across ROIs ──────────────────────────────
    src_r2_list <- list()
    for (this_roi in all_rois) {
        src_fitted_path <- file.path(source_out_dir, this_roi,
                                     paste0("trial_level_fitted_", this_roi, ".csv"))
        src_fitted <- read_safe(src_fitted_path)
        if (is.null(src_fitted)) next
        fit_col_src <- paste0(src_model, "_fitted")
        if (!fit_col_src %in% names(src_fitted)) next
        if (!"pain_rating" %in% names(src_fitted)) next
        
        src_r2_list[[this_roi]] <- per_subject_r2(src_fitted, "pain_rating", fit_col_src) %>%
            mutate(roi = this_roi)
    }
    
    if (length(src_r2_list) == 0L) {
        message("  [SKIP] No source fitted values found for: ", src_model)
        next
    }
    
    # Use the ROI with the highest mean R² as the source comparator
    src_r2_means <- map_dbl(src_r2_list, ~mean(.x$r2, na.rm = TRUE))
    best_roi      <- names(which.max(src_r2_means))
    src_r2        <- src_r2_list[[best_roi]] %>%
        rename(r2_source = r2) %>%
        select(subjid_uid, r2_source)
    
    joined <- inner_join(ch_r2, src_r2, by = "subjid_uid")
    
    tt_row <- run_paired_ttest(joined$r2_source, joined$r2_channel, concept) %>%
        mutate(concept = concept, channel_model = ch_model,
               source_model = src_model, best_roi = best_roi,
               mean_r2_channel = mean(joined$r2_channel, na.rm = TRUE),
               mean_r2_source  = mean(joined$r2_source,  na.rm = TRUE))
    
    test2_rows[[i]] <- tt_row
    
    test2_plot_data[[i]] <- joined %>%
        pivot_longer(c(r2_channel, r2_source),
                     names_to = "domain", values_to = "r2") %>%
        mutate(concept   = concept,
               domain    = recode(domain,
                                  "r2_channel" = "Channel",
                                  "r2_source"  = paste0("Source (", best_roi, ")")))
}

test2_results <- bind_rows(test2_rows) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr"),
           sig   = sig_label(p_adj))

write_csv(test2_results, file.path(out_dir, "test2_fit_ttest.csv"))
message("  Saved test2_fit_ttest.csv  (", nrow(test2_results), " rows)")

# ── Test 2 plot: paired dot-and-line ─────────────────────────────────────────
test2_plot_df <- bind_rows(test2_plot_data)

if (nrow(test2_plot_df) > 0) {
    # Build explicit paired data: one row per (subject, concept, domain)
    # Pivot wide then long to get two-column structure for geom_line grouping
    test2_paired <- test2_plot_df %>%
        pivot_wider(names_from = domain, values_from = r2,
                    id_cols = c(subjid_uid, concept)) %>%
        pivot_longer(
            cols      = -c(subjid_uid, concept),
            names_to  = "domain",
            values_to = "r2"
        ) %>%
        mutate(domain = factor(domain))   # consistent ordering
    
    domain_lvls <- levels(test2_paired$domain)
    
    p2 <- ggplot(test2_paired,
                 aes(x = domain, y = r2, colour = domain)) +
        geom_line(aes(group = subjid_uid),
                  colour = "grey60", alpha = 0.4, linewidth = 0.35) +
        geom_point(size = 1.4, alpha = 0.7) +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 4,
                     colour = "black") +
        stat_summary(fun.data = mean_se, geom = "errorbar",
                     width = 0.15, colour = "black") +
        stat_compare_means(method = "t.test", paired = TRUE,
                           label = "p.signif", size = 3.5,
                           comparisons = list(domain_lvls[1:2])) +
        facet_wrap(~concept, scales = "free_y", nrow = 1) +
        labs(title    = "Test 2: Per-Subject R² — Channel vs Source GAMM",
             subtitle = "Lines connect the same subject; diamond = mean ± SE",
             x = NULL, y = "Per-subject R²", colour = NULL) +
        theme_minimal(base_size = 10) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 20, hjust = 1))
    
    ggsave(file.path(out_dir, "test2_fit_plot.png"),
           p2, width = max(8, nrow(concept_pairs) * 3), height = 5, dpi = 200)
    message("  Saved test2_fit_plot.png")
}

# =============================================================================
# TEST 3: REPEATED-MEASURES ANOVA  —  METRICS ACROSS ROIs
# =============================================================================
message("\n", strrep("=", 60))
message("TEST 3: rmANOVA — BI_pre / CoG_pre / delta_ERD across ROIs")
message(strrep("=", 60))

test3_plot_list <- list()

for (metric in rmanova_metrics) {
    if (!metric %in% names(src_subject_ga)) {
        message("  [SKIP] '", metric, "' not in subject GA data.")
        next
    }
    message("  Running rmANOVA for: ", metric)
    
    # Need at least 2 ROIs and 3 subjects with complete data across all ROIs
    rm_df <- src_subject_ga %>%
        select(subjid_uid, roi, all_of(metric)) %>%
        rename(value = all_of(metric)) %>%
        filter(!is.na(value)) %>%
        # Keep only subjects with data in all ROIs
        group_by(subjid_uid) %>%
        filter(n_distinct(roi) == length(all_rois)) %>%
        ungroup() %>%
        mutate(
            subjid_uid = as.factor(subjid_uid),
            roi        = as.factor(roi)
        )
    
    if (n_distinct(rm_df$subjid_uid) < 3L) {
        message("  [SKIP] Fewer than 3 complete subjects for: ", metric)
        next
    }
    
    # ── Fit rmANOVA ─────────────────────────────────────────────────────────
    anova_tbl <- tryCatch({
        if (.afex_ok) {
            fit <- afex::aov_ez(id = "subjid_uid", dv = "value",
                                data = rm_df, within = "roi",
                                anova_table = list(correction = "GG"))
            as_tibble(fit$anova_table, rownames = "Effect") %>%
                rename(
                    num_df   = `num Df`,
                    den_df   = `den Df`,
                    MSE      = MSE,
                    F_stat   = F,
                    p_value  = `Pr(>F)`,
                    p_adj_GG = `Pr(>F)`   # GG-corrected p is reported directly by afex
                ) %>%
                mutate(metric = metric, correction = "Greenhouse-Geisser")
        } else {
            # Fallback: base R aov with Error()
            fit  <- aov(value ~ roi + Error(subjid_uid / roi), data = rm_df)
            s    <- summary(fit)$`Error: subjid_uid:roi`[[1]]
            tibble(
                Effect    = "roi",
                num_df    = s["roi", "Df"],
                den_df    = NA,
                MSE       = NA,
                F_stat    = s["roi", "F value"],
                p_value   = s["roi", "Pr(>F)"],
                metric    = metric,
                correction = "Sphericity assumed (install afex for GG)"
            )
        }
    }, error = function(e) {
        message("  [WARN] rmANOVA failed for ", metric, ": ", conditionMessage(e))
        NULL
    })
    
    if (is.null(anova_tbl)) next
    
    anova_tbl$sig <- sig_label(anova_tbl$p_value)
    
    # ── Post-hoc pairwise (FDR corrected) ───────────────────────────────────
    posthoc <- pairwise_fdr(rm_df, "value", "roi", "subjid_uid")
    
    out_df <- list(anova = anova_tbl, posthoc = posthoc)
    write_csv(anova_tbl, file.path(out_dir, paste0("test3_rmanova_", metric, ".csv")))
    write_csv(posthoc,   file.path(out_dir, paste0("test3_posthoc_", metric, ".csv")))
    message("  Saved test3_rmanova_", metric, ".csv + test3_posthoc_", metric, ".csv")
    
    test3_plot_list[[metric]] <- rm_df %>% mutate(metric = metric)
}

# ── Test 3 plot: faceted boxplot per metric × ROI ────────────────────────────
if (length(test3_plot_list) > 0L) {
    test3_plot_df <- bind_rows(test3_plot_list)
    
    p3 <- ggplot(test3_plot_df, aes(x = roi, y = value, fill = roi)) +
        geom_boxplot(outlier.size = 0.5, alpha = 0.8, width = 0.6) +
        geom_jitter(width = 0.12, size = 0.8, alpha = 0.4) +
        facet_wrap(~metric, scales = "free_y", ncol = length(rmanova_metrics)) +
        labs(title  = "Test 3: Source Metrics Across ROIs (rmANOVA)",
             subtitle = "Each point = one subject GA mean; boxes = group distribution",
             x = "ROI", y = "Metric value", fill = "ROI") +
        theme_minimal(base_size = 10) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 30, hjust = 1))
    
    ggsave(file.path(out_dir, "test3_rmanova_plot.png"),
           p3, width = length(rmanova_metrics) * 4 + 2, height = 5, dpi = 200)
    message("  Saved test3_rmanova_plot.png")
}

# =============================================================================
# TEST 4: ONE-WAY ANOVA  —  SOURCE GAMM DEVIANCE EXPLAINED ACROSS ROIs
# =============================================================================
message("\n", strrep("=", 60))
message("TEST 4: ANOVA — source GAMM deviance across ROIs")
message(strrep("=", 60))

test4_rows  <- list()
test4_plot_list <- list()

for (src_model in source_models_for_r2) {
    r2_by_roi <- list()
    
    for (this_roi in all_rois) {
        src_fitted_path <- file.path(source_out_dir, this_roi,
                                     paste0("trial_level_fitted_", this_roi, ".csv"))
        src_fitted <- read_safe(src_fitted_path)
        if (is.null(src_fitted)) next
        fit_col <- paste0(src_model, "_fitted")
        if (!fit_col %in% names(src_fitted) || !"pain_rating" %in% names(src_fitted)) next
        
        r2_sub <- per_subject_r2(src_fitted, "pain_rating", fit_col) %>%
            mutate(roi = this_roi, model = src_model)
        r2_by_roi[[this_roi]] <- r2_sub
    }
    
    if (length(r2_by_roi) < 2L) {
        message("  [SKIP] Fewer than 2 ROIs for model: ", src_model)
        next
    }
    
    anova_df <- bind_rows(r2_by_roi) %>%
        filter(!is.na(r2)) %>%
        mutate(roi = as.factor(roi))
    
    if (n_distinct(anova_df$roi) < 2L || nrow(anova_df) < 6L) next
    
    anova_fit  <- aov(r2 ~ roi, data = anova_df)
    anova_summ <- broom::tidy(anova_fit) %>%
        mutate(model = src_model, sig = sig_label(p.value))
    
    posthoc4 <- pairwise_fdr(anova_df, "r2", "roi", "subjid_uid")
    
    test4_rows[[src_model]] <- anova_summ
    write_csv(anova_summ, file.path(out_dir, paste0("test4_roi_anova_", src_model, ".csv")))
    write_csv(posthoc4,   file.path(out_dir, paste0("test4_posthoc_", src_model, ".csv")))
    
    test4_plot_list[[src_model]] <- anova_df
    message("  Saved test4_roi_anova_", src_model, ".csv")
}

# ── Test 4 plot: mean R² per ROI per model ────────────────────────────────────
if (length(test4_plot_list) > 0L) {
    test4_plot_df <- bind_rows(test4_plot_list)
    
    test4_summary <- test4_plot_df %>%
        group_by(model, roi) %>%
        summarise(mean_r2 = mean(r2, na.rm = TRUE),
                  se_r2   = sd(r2, na.rm = TRUE) / sqrt(sum(!is.na(r2))),
                  .groups = "drop")
    
    p4 <- ggplot(test4_summary, aes(x = roi, y = mean_r2, fill = roi)) +
        geom_col(width = 0.65, alpha = 0.85) +
        geom_errorbar(aes(ymin = mean_r2 - se_r2, ymax = mean_r2 + se_r2),
                      width = 0.2) +
        facet_wrap(~model, scales = "free_y") +
        labs(title   = "Test 4: Per-Subject R² Across ROIs by Model",
             subtitle = "Mean ± SE across subjects; one-way ANOVA tests ROI differences",
             x = "ROI", y = "Mean per-subject R²", fill = "ROI") +
        theme_minimal(base_size = 10) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 30, hjust = 1))
    
    ggsave(file.path(out_dir, "test4_roi_anova_plot.png"),
           p4, width = length(source_models_for_r2) * 3 + 2, height = 5, dpi = 200)
    message("  Saved test4_roi_anova_plot.png")
}

# =============================================================================
# TEST 5: RAYLEIGH TEST  —  PHASE UNIFORMITY
# =============================================================================
# Self-contained implementation; no circular package required.
#
# Rayleigh statistic:  z = n * R²
#   R = mean resultant length = |Σ exp(iθ)| / n
# P-value approximation: p ≈ exp(-z) * (1 + (2z - z²)/(4n) - (24z - 132z² + 76z³ - 9z⁴)/(288n²))
# This is the Mardia-Jupp approximation, accurate for n ≥ 10.
# =============================================================================
message("\n", strrep("=", 60))
message("TEST 5: Rayleigh test — phase uniformity")
message(strrep("=", 60))

rayleigh_test <- function(phases) {
    phases <- phases[!is.na(phases)]
    n      <- length(phases)
    if (n < 5L) return(list(n = n, R = NA, z = NA, p = NA))
    
    C   <- mean(cos(phases))
    S   <- mean(sin(phases))
    R   <- sqrt(C^2 + S^2)       # mean resultant length ∈ [0, 1]
    z   <- n * R^2
    # Mardia-Jupp p-value approximation
    p   <- exp(-z) * (1 + (2*z - z^2) / (4*n) -
                          (24*z - 132*z^2 + 76*z^3 - 9*z^4) / (288*n^2))
    p   <- max(0, min(1, p))     # clamp to [0, 1]
    list(n = n, R = round(R, 4), z = round(z, 3), p = p,
         mean_phase = atan2(S, C))
}

phase_cols <- intersect(c("slow_phase", "slow_phase_post"), names(src_master))

for (pcol in phase_cols) {
    message("  Phase column: ", pcol)
    ray_rows <- list()
    
    for (this_roi in all_rois) {
        roi_trials <- src_master %>%
            filter(roi == this_roi, !is.na(.data[[pcol]])) %>%
            select(subjid_uid, all_of(pcol))
        
        for (subj in unique(roi_trials$subjid_uid)) {
            phases <- roi_trials %>% filter(subjid_uid == subj) %>% pull(all_of(pcol))
            rt     <- rayleigh_test(phases)
            ray_rows[[length(ray_rows) + 1L]] <- tibble(
                subjid_uid  = subj,
                roi         = this_roi,
                phase_col   = pcol,
                n_trials    = rt$n,
                R_mean      = rt$R,
                rayleigh_z  = rt$z,
                p_value     = rt$p,
                mean_phase  = rt$mean_phase %||% NA
            )
        }
    }
    
    ray_df <- bind_rows(ray_rows) %>%
        mutate(
            p_adj = p.adjust(p_value, method = "fdr"),
            sig   = sig_label(p_adj)
        )
    
    write_csv(ray_df, file.path(out_dir, paste0("test5_rayleigh_", pcol, ".csv")))
    message("  Saved test5_rayleigh_", pcol, ".csv  (", nrow(ray_df), " rows)")
    
    # ROI-level summary
    ray_summary <- ray_df %>%
        group_by(roi) %>%
        summarise(
            n_subjects    = n(),
            n_sig_fdr     = sum(p_adj < alpha_level, na.rm = TRUE),
            pct_sig       = round(100 * n_sig_fdr / n_subjects, 1),
            mean_R        = round(mean(R_mean,    na.rm = TRUE), 3),
            median_R      = round(median(R_mean,  na.rm = TRUE), 3),
            mean_rayleigh_z = round(mean(rayleigh_z, na.rm = TRUE), 2),
            .groups = "drop"
        ) %>%
        mutate(phase_col = pcol)
    
    write_csv(ray_summary,
              file.path(out_dir, paste0("test5_rayleigh_summary_", pcol, ".csv")))
    
    # ── Test 5 plot: mean resultant length heatmap + significance overlay ──────
    heatmap_df <- ray_df %>%
        select(subjid_uid, roi, R_mean, sig) %>%
        mutate(roi = factor(roi, levels = all_rois))
    
    n_subj    <- n_distinct(heatmap_df$subjid_uid)
    subj_order <- heatmap_df %>%
        group_by(subjid_uid) %>%
        summarise(mean_R = mean(R_mean, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_R)) %>%
        pull(subjid_uid)
    
    heatmap_df <- heatmap_df %>%
        mutate(subjid_uid = factor(subjid_uid, levels = subj_order))
    
    p5 <- ggplot(heatmap_df, aes(x = roi, y = subjid_uid, fill = R_mean)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        geom_text(aes(label = sig), size = 2.5, colour = "black") +
        scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1),
                            name = "R\n(mean resultant)") +
        labs(
            title    = paste0("Test 5: Rayleigh Test — ", pcol),
            subtitle = "Colour = mean resultant length R; label = FDR significance",
            x = "ROI", y = "Subject"
        ) +
        theme_minimal(base_size = 9) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              axis.text.y = element_text(size = 6))
    
    ggsave(file.path(out_dir, paste0("test5_rayleigh_plot_", pcol, ".png")),
           p5, width = max(6, length(all_rois) * 1.5 + 2),
           height = max(5, n_subj * 0.3 + 2), dpi = 200)
    message("  Saved test5_rayleigh_plot_", pcol, ".png")
}

# =============================================================================
# SUMMARY REPORT
# =============================================================================
message("\nWriting summary report...")

report <- c(
    "CLASSICAL TESTS SUMMARY REPORT",
    strrep("=", 60),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
    ""
)

# Test 1
report <- c(report, "TEST 1: ERD Significance (pre vs post alpha power)",
            strrep("-", 40))
if (exists("test1_results") && nrow(test1_results) > 0) {
    for (i in seq_len(nrow(test1_results))) {
        r <- test1_results[i, ]
        report <- c(report, sprintf(
            "  %-35s  t(%s)=%.2f  p(adj)=%.4f  d=%.2f  %s",
            r$label, r$df, r$t %||% NA, r$p_adj %||% NA,
            r$cohens_d %||% NA, r$sig
        ))
    }
} else {
    report <- c(report, "  No results produced.")
}

# Test 2
report <- c(report, "", "TEST 2: Channel vs Source R² (paired t-test)",
            strrep("-", 40))
if (exists("test2_results") && nrow(test2_results) > 0) {
    for (i in seq_len(nrow(test2_results))) {
        r <- test2_results[i, ]
        winner <- if (!is.na(r$mean_diff) && r$mean_diff > 0) "SOURCE" else "CHANNEL"
        report <- c(report, sprintf(
            "  %-20s  t(%s)=%.2f  p(adj)=%.4f  %s  (mean R² chan=%.3f src=%.3f)",
            r$concept, r$df, r$t %||% NA, r$p_adj %||% NA,
            r$sig, r$mean_r2_channel %||% NA, r$mean_r2_source %||% NA
        ))
    }
}

# Test 3
report <- c(report, "", "TEST 3: rmANOVA — metrics across ROIs",
            strrep("-", 40))
for (metric in rmanova_metrics) {
    f <- file.path(out_dir, paste0("test3_rmanova_", metric, ".csv"))
    if (file.exists(f)) {
        tbl <- read_csv(f, show_col_types = FALSE)
        row <- tbl[1, ]
        report <- c(report, sprintf(
            "  %-15s  F(%s,%s)=%.2f  p=%.4f  %s",
            metric,
            round(row$num_df %||% NA, 1),
            round(row$den_df %||% NA, 1),
            row$F_stat %||% NA,
            row$p_value %||% NA,
            row$sig %||% ""
        ))
    }
}

# Test 4
report <- c(report, "", "TEST 4: ANOVA — deviance across ROIs",
            strrep("-", 40))
for (src_model in source_models_for_r2) {
    f <- file.path(out_dir, paste0("test4_roi_anova_", src_model, ".csv"))
    if (file.exists(f)) {
        tbl <- read_csv(f, show_col_types = FALSE) %>% filter(term == "roi")
        res <- read_csv(f, show_col_types = FALSE) %>% filter(term == "Residuals")
        if (nrow(tbl) > 0) {
            report <- c(report, sprintf(
                "  %-25s  F(%s,%s)=%.2f  p=%.4f  %s",
                src_model,
                round(tbl$df[1], 1),
                if (nrow(res) > 0) round(res$df[1], 1) else "?",
                tbl$statistic[1] %||% NA,
                tbl$p.value[1]   %||% NA,
                sig_label(tbl$p.value[1])
            ))
        }
    }
}

# Test 5
report <- c(report, "", "TEST 5: Rayleigh test — phase uniformity",
            strrep("-", 40))
for (pcol in phase_cols) {
    f <- file.path(out_dir, paste0("test5_rayleigh_summary_", pcol, ".csv"))
    if (file.exists(f)) {
        smry <- read_csv(f, show_col_types = FALSE)
        report <- c(report, paste0("  ", pcol, ":"))
        for (i in seq_len(nrow(smry))) {
            r <- smry[i, ]
            report <- c(report, sprintf(
                "    %-15s  mean_R=%.3f  n_sig=%d/%d (%.1f%%)",
                r$roi, r$mean_R, r$n_sig_fdr, r$n_subjects, r$pct_sig
            ))
        }
    }
}

writeLines(report, file.path(out_dir, "classical_tests_report.txt"))
message("\nClassical tests complete. Outputs: ", out_dir)
message(paste(report, collapse = "\n"))