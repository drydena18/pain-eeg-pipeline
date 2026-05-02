# =============================================================================
# compare_channel_source_gamm.R
# -----------------------------------------------------------------------------
# Systematically compares channel-wise and source-wise GAMMs across three
# dimensions:
#
#   1. Model fit quality    AIC / BIC / deviance explained / R²
#   2. Effect significance  p-values and EDF (effective degrees of freedom)
#                           for matched smooth terms
#   3. Effect concordance   Direction and shape agreement between channel
#                           and source smooth estimates for the same metric
#
# Inputs (read automatically; no manual edits needed)
# ────────────────────────────────────────────────────
#   <channel_out_dir>/model_comparison_v2.csv         channel fit table
#   <channel_out_dir>/<model>.rds                     fitted channel models
#   <channel_out_dir>/trial_level_fitted_values_v2.csv
#
#   <source_out_dir>/model_comparison_all_rois.csv    source fit table
#   <source_out_dir>/<roi>/<model>.rds                fitted source models (per ROI)
#   <source_out_dir>/<roi>/trial_level_fitted_<roi>.csv
#
# Outputs  (<compare_out_dir>/)
# ──────────────────────────────
#   fit_quality_comparison.csv        AIC/deviance side-by-side per matched model
#   effect_significance_comparison.csv  term-level p + EDF for channel vs source
#   smooth_concordance.csv            correlation of partial effects per metric
#   fit_quality_plot.png              bar chart — deviance explained by domain
#   delta_aic_plot.png                channel vs best source ROI Δ AIC per model
#   smooth_overlay_<metric>.png       channel vs source smooth estimate overlay
#   concordance_heatmap.png           ROI × metric concordance (r) heatmap
#   comparison_report.txt             plain-text narrative summary
# =============================================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(mgcv)
library(ggplot2)
library(tibble)

# %||% is native in R >= 4.4; define a fallback for older installations
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[[1]])) a else b
}

# =============================================================================
# USER SETTINGS
# =============================================================================
channel_out_dir <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs_v2"
source_out_dir  <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs_source"
compare_out_dir <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_comparison"

dir.create(compare_out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# CONCEPT MAP
# Channel model name -> (concept label, source model name, predictor term)
# Used to align models across domains for fair comparison.
# =============================================================================
concept_map <- tribble(
  ~channel_model,                  ~concept,           ~source_model,       ~term_z,
  "m00_baseline",                  "Baseline",          "m00_baseline",      NA,
  "m02_sf_logratio",               "Log ratio",         "m02_lr_pre",        "LR_pre_z",
  "m03_sf_balance",                "Balance index",     "m01_bi_pre",        "BI_pre_z",
  "m04_slow_frac",                 "Slow fraction",     "m03_cog_pre",       "CoG_pre_z",
  "m05_tensor_slow_fast",          "Tensor slow×fast",  "m05_tensor_slow_fast", NA,
  "m06_sf_logratio_aperiodic",     "LR (FOOOF ctrl)",   "m14_fooof_lr_pre",  "LR_pre_z",
  "m07_sf_balance_aperiodic",      "BI (FOOOF ctrl)",   "m13_fooof_bi_pre",  "BI_pre_z"
)

# =============================================================================
# HELPERS
# =============================================================================
load_rds_safe <- function(path) {
  tryCatch(readRDS(path), error = function(e) {
    message("[WARN] Could not load RDS: ", path, " — ", conditionMessage(e))
    NULL
  })
}

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read_csv(path, show_col_types = FALSE),
           error = function(e) { message("[WARN] ", path); NULL })
}

# Extract smooth term summary from a fitted bam/gam.
# Returns a tibble with one row per smooth term.
extract_smooth_summary <- function(mod, model_nm, domain, roi = NA_character_) {
  if (is.null(mod)) return(NULL)
  sm  <- summary(mod)
  if (is.null(sm$s.table)) return(NULL)
  as_tibble(sm$s.table, rownames = "term") %>%
    rename(edf = edf, ref_df = `Ref.df`, F_stat = F, p_value = `p-value`) %>%
    mutate(model_name = model_nm, domain = domain, roi = roi)
}

# Get the partial effect of a smooth over its range (for overlay plots).
# Returns a tibble: x (predictor value), fit, se.
get_smooth_effect <- function(mod, term_z) {
  if (is.null(mod) || is.null(term_z) || is.na(term_z)) return(NULL)
  tryCatch({
    # Use plot.gam internals via gratia-style approach
    nd_range <- range(mod$model[[term_z]], na.rm = TRUE)
    nd <- data.frame(x = seq(nd_range[1], nd_range[2], length.out = 200))
    names(nd) <- term_z
    
    # Fill all other predictors at their reference level
    for (v in names(mod$model)) {
      if (v == term_z || v == "pain_rating") next
      if (is.factor(mod$model[[v]])) {
        nd[[v]] <- factor(levels(mod$model[[v]])[1], levels = levels(mod$model[[v]]))
      } else {
        nd[[v]] <- 0  # z-scored predictors: 0 = mean
      }
    }
    
    pred  <- predict(mod, newdata = nd, type = "terms", se.fit = TRUE,
                     terms = paste0("s(", term_z, ")"))
    tibble(
      x_z   = nd[[term_z]],
      fit   = as.numeric(pred$fit),
      se    = as.numeric(pred$se.fit),
      lo    = fit - 1.96 * se,
      hi    = fit + 1.96 * se
    )
  }, error = function(e) NULL)
}

# =============================================================================
# 1. LOAD FIT QUALITY TABLES
# =============================================================================
message("Loading fit quality tables...")

channel_comp <- read_csv_safe(file.path(channel_out_dir, "model_comparison_v2.csv"))
source_comp  <- read_csv_safe(file.path(source_out_dir,  "model_comparison_all_rois.csv"))

if (is.null(channel_comp)) stop("Channel comparison table not found.")
if (is.null(source_comp))  stop("Source comparison table not found.")

all_rois <- sort(unique(source_comp$roi))
message("Source ROIs: ", paste(all_rois, collapse = ", "))

# =============================================================================
# 2. FIT QUALITY COMPARISON TABLE
# =============================================================================
message("Building fit quality comparison...")

# For each concept, find the best source ROI (lowest AIC for that model)
source_best <- source_comp %>%
  inner_join(concept_map %>% select(concept, source_model),
             by = c("model_name" = "source_model")) %>%
  group_by(concept, model_name) %>%
  slice_min(AIC, n = 1L) %>%
  ungroup() %>%
  rename(
    source_model   = model_name,
    source_AIC     = AIC,
    source_BIC     = BIC,
    source_devexpl = dev_expl,
    source_r_sq    = r_sq,
    source_n       = n,
    best_roi       = roi
  )

channel_matched <- channel_comp %>%
  inner_join(concept_map %>% select(concept, channel_model),
             by = c("model_name" = "channel_model")) %>%
  rename(
    channel_model   = model_name,
    channel_AIC     = AIC,
    channel_BIC     = BIC,
    channel_devexpl = dev_expl,
    channel_n       = n
  ) %>%
  select(concept, channel_model, channel_AIC, channel_BIC,
         channel_devexpl, channel_n)

fit_quality <- channel_matched %>%
  left_join(source_best, by = "concept") %>%
  mutate(
    delta_AIC_chan_minus_src = channel_AIC - source_AIC,
    delta_devexpl            = source_devexpl - channel_devexpl
  ) %>%
  arrange(concept)

write_csv(fit_quality, file.path(compare_out_dir, "fit_quality_comparison.csv"))
message("  Saved fit_quality_comparison.csv")

# =============================================================================
# 3. EFFECT SIGNIFICANCE COMPARISON
# =============================================================================
message("Extracting smooth-term summaries from fitted models...")

# Channel models
channel_smooth_rows <- list()
for (i in seq_len(nrow(concept_map))) {
  cm  <- concept_map$channel_model[i]
  rds <- file.path(channel_out_dir, paste0(cm, ".rds"))
  mod <- load_rds_safe(rds)
  channel_smooth_rows[[i]] <- extract_smooth_summary(mod, cm, "channel", NA)
}
channel_smooth <- bind_rows(channel_smooth_rows)

# Source models (iterate over ROIs × concepts)
source_smooth_rows <- list()
idx <- 1L
for (this_roi in all_rois) {
  for (i in seq_len(nrow(concept_map))) {
    sm  <- concept_map$source_model[i]
    rds <- file.path(source_out_dir, this_roi, paste0(sm, ".rds"))
    mod <- load_rds_safe(rds)
    source_smooth_rows[[idx]] <- extract_smooth_summary(mod, sm, "source", this_roi)
    idx <- idx + 1L
  }
}
source_smooth <- bind_rows(source_smooth_rows)

all_smooth <- bind_rows(channel_smooth, source_smooth)
write_csv(all_smooth, file.path(compare_out_dir, "effect_significance_comparison.csv"))
message("  Saved effect_significance_comparison.csv")

# =============================================================================
# 4. SMOOTH CONCORDANCE  (channel vs source partial effect correlation)
# =============================================================================
message("Computing smooth concordance...")

concordance_rows <- list()
conc_idx <- 1L

for (i in seq_len(nrow(concept_map))) {
  term_z  <- concept_map$term_z[i]
  concept <- concept_map$concept[i]
  cm      <- concept_map$channel_model[i]
  sm      <- concept_map$source_model[i]
  
  if (is.na(term_z)) next
  
  ch_mod <- load_rds_safe(file.path(channel_out_dir, paste0(cm, ".rds")))
  ch_eff <- get_smooth_effect(ch_mod, term_z)
  if (is.null(ch_eff)) next
  
  for (this_roi in all_rois) {
    src_mod <- load_rds_safe(
      file.path(source_out_dir, this_roi, paste0(sm, ".rds"))
    )
    src_eff <- get_smooth_effect(src_mod, term_z)
    if (is.null(src_eff)) next
    
    # Interpolate source onto channel x-grid for correlation
    src_fit_interp <- approx(src_eff$x_z, src_eff$fit, xout = ch_eff$x_z,
                             rule = 2)$y
    r <- cor(ch_eff$fit, src_fit_interp, use = "complete.obs")
    
    concordance_rows[[conc_idx]] <- tibble(
      concept = concept, term_z = term_z, roi = this_roi,
      r_concordance = r,
      same_direction = sign(r) == 1L
    )
    conc_idx <- conc_idx + 1L
  }
}

concordance <- bind_rows(concordance_rows)
write_csv(concordance, file.path(compare_out_dir, "smooth_concordance.csv"))
message("  Saved smooth_concordance.csv")

# =============================================================================
# 5. PLOTS
# =============================================================================
message("Generating comparison plots...")

# ── 5a. Deviance explained bar chart ─────────────────────────────────────────
if (nrow(fit_quality) > 0) {
  fq_long <- fit_quality %>%
    select(concept, channel_devexpl, source_devexpl, best_roi) %>%
    pivot_longer(c(channel_devexpl, source_devexpl),
                 names_to = "domain", values_to = "dev_expl") %>%
    mutate(
      domain = recode(domain,
                      "channel_devexpl" = "Channel",
                      "source_devexpl"  = "Source"),
      label = if_else(domain == "Source",
                      paste0("Source\n(", best_roi, ")"), domain)
    )
  
  p_dev <- ggplot(fq_long, aes(x = reorder(concept, dev_expl), y = dev_expl,
                               fill = domain)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = c("Channel" = "steelblue", "Source" = "tomato")) +
    labs(title = "Deviance Explained: Channel vs Source GAMM",
         x = NULL, y = "Deviance explained", fill = "Domain") +
    theme_minimal(base_size = 11)
  
  ggsave(file.path(compare_out_dir, "fit_quality_plot.png"),
         p_dev, width = 8, height = 5, dpi = 200)
  
  # ── 5b. Δ AIC plot (positive = source beats channel) ─────────────────────
  p_aic <- ggplot(fit_quality,
                  aes(x = reorder(concept, delta_AIC_chan_minus_src),
                      y = delta_AIC_chan_minus_src,
                      fill = delta_AIC_chan_minus_src > 0)) +
    geom_col(width = 0.65) +
    coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue"),
                      labels = c("TRUE"  = "Source better",
                                 "FALSE" = "Channel better")) +
    labs(title = "Δ AIC (Channel − Source): positive = source fits better",
         x = NULL, y = "Δ AIC (channel − source)", fill = NULL) +
    theme_minimal(base_size = 11)
  
  ggsave(file.path(compare_out_dir, "delta_aic_plot.png"),
         p_aic, width = 8, height = 5, dpi = 200)
}

# ── 5c. Smooth overlay: channel vs all source ROIs per metric ─────────────
unique_terms <- concept_map %>% filter(!is.na(term_z)) %>%
  distinct(concept, channel_model, source_model, term_z)

for (i in seq_len(nrow(unique_terms))) {
  term_z  <- unique_terms$term_z[i]
  concept <- unique_terms$concept[i]
  cm      <- unique_terms$channel_model[i]
  sm      <- unique_terms$source_model[i]
  
  ch_mod <- load_rds_safe(file.path(channel_out_dir, paste0(cm, ".rds")))
  ch_eff <- get_smooth_effect(ch_mod, term_z)
  
  eff_list <- list()
  if (!is.null(ch_eff)) {
    eff_list[["Channel"]] <- ch_eff %>% mutate(roi = "Channel")
  }
  
  for (this_roi in all_rois) {
    src_mod <- load_rds_safe(
      file.path(source_out_dir, this_roi, paste0(sm, ".rds"))
    )
    src_eff <- get_smooth_effect(src_mod, term_z)
    if (!is.null(src_eff)) {
      eff_list[[this_roi]] <- src_eff %>% mutate(roi = this_roi)
    }
  }
  
  if (length(eff_list) < 2L) next
  
  eff_df <- bind_rows(eff_list) %>%
    mutate(is_channel = roi == "Channel")
  
  p_smooth <- ggplot(eff_df, aes(x = x_z, y = fit,
                                 colour = roi, linetype = is_channel,
                                 fill = roi)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4, linetype = "dashed") +
    scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                          guide = "none") +
    labs(title = paste0(concept, "  —  smooth partial effect (", term_z, ")"),
         subtitle = "Dashed = channel; solid = source ROIs",
         x = paste0(term_z, "  (z-scored)"),
         y = "Partial effect on pain rating",
         colour = "Domain / ROI", fill = "Domain / ROI") +
    theme_minimal(base_size = 11)
  
  out_name <- paste0("smooth_overlay_", str_replace_all(concept, "\\s+", "_"), ".png")
  ggsave(file.path(compare_out_dir, out_name),
         p_smooth, width = 8, height = 5, dpi = 200)
}

# ── 5d. Concordance heatmap (ROI × metric) ───────────────────────────────
if (nrow(concordance) > 0) {
  p_heatmap <- ggplot(concordance,
                      aes(x = roi, y = concept, fill = r_concordance)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", r_concordance)),
              size = 2.8, colour = "black") +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                         midpoint = 0, limits = c(-1, 1),
                         name = "r (concordance)") +
    labs(title = "Smooth Effect Concordance: Channel vs Source ROI",
         subtitle = "r = correlation of partial effect curves",
         x = "Source ROI", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(file.path(compare_out_dir, "concordance_heatmap.png"),
         p_heatmap, width = max(6, length(all_rois) * 1.2 + 2), height = 6,
         dpi = 200)
}

# =============================================================================
# 6. PLAIN-TEXT SUMMARY REPORT
# =============================================================================
message("Writing summary report...")

report_lines <- c(
  "CHANNEL vs SOURCE GAMM COMPARISON REPORT",
  strrep("=", 60),
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "1. FIT QUALITY (deviance explained)",
  strrep("-", 40)
)

if (nrow(fit_quality) > 0) {
  for (i in seq_len(nrow(fit_quality))) {
    r <- fit_quality[i, ]
    winner <- if (!is.na(r$delta_AIC_chan_minus_src) &&
                  r$delta_AIC_chan_minus_src > 2) "SOURCE" else
                    if (!is.na(r$delta_AIC_chan_minus_src) &&
                        r$delta_AIC_chan_minus_src < -2) "CHANNEL" else "TIED"
    report_lines <- c(report_lines, sprintf(
      "  %-25s  Chan devexpl=%.3f  Src devexpl=%.3f  ΔAIC=%.1f  → %s (best ROI: %s)",
      r$concept,
      r$channel_devexpl %||% NA,
      r$source_devexpl  %||% NA,
      r$delta_AIC_chan_minus_src %||% NA,
      winner,
      r$best_roi %||% "?"
    ))
  }
}

report_lines <- c(report_lines, "",
                  "2. SMOOTH EFFECT CONCORDANCE",
                  strrep("-", 40)
)
if (nrow(concordance) > 0) {
  conc_summary <- concordance %>%
    group_by(concept, term_z) %>%
    summarise(
      mean_r    = mean(r_concordance, na.rm = TRUE),
      n_agree   = sum(same_direction, na.rm = TRUE),
      n_roi     = n(),
      .groups   = "drop"
    )
  for (i in seq_len(nrow(conc_summary))) {
    r <- conc_summary[i, ]
    report_lines <- c(report_lines, sprintf(
      "  %-20s  mean r=%.2f  direction agreement: %d / %d ROIs",
      r$concept, r$mean_r, r$n_agree, r$n_roi
    ))
  }
}

report_lines <- c(report_lines, "",
                  "3. DOMAIN RECOMMENDATION",
                  strrep("-", 40)
)

if (nrow(fit_quality) > 0) {
  n_src_better <- sum(fit_quality$delta_AIC_chan_minus_src > 2, na.rm = TRUE)
  n_ch_better  <- sum(fit_quality$delta_AIC_chan_minus_src < -2, na.rm = TRUE)
  n_tied       <- nrow(fit_quality) - n_src_better - n_ch_better
  
  report_lines <- c(report_lines,
                    sprintf("  Source better (ΔAIC > 2) : %d / %d models", n_src_better, nrow(fit_quality)),
                    sprintf("  Channel better (ΔAIC < -2): %d / %d models", n_ch_better, nrow(fit_quality)),
                    sprintf("  Tied (|ΔAIC| ≤ 2)         : %d / %d models", n_tied, nrow(fit_quality))
  )
  if (n_src_better > n_ch_better) {
    report_lines <- c(report_lines, "",
                      "  Overall: SOURCE-SPACE metrics provide better GAMM fit.",
                      "  Interpretation: localised source reconstruction reduces",
                      "  noise and improves the alpha-pain relationship signal.")
  } else if (n_ch_better > n_src_better) {
    report_lines <- c(report_lines, "",
                      "  Overall: CHANNEL-SPACE metrics provide better GAMM fit.",
                      "  Interpretation: the spatial mixing at the scalp level may",
                      "  capture integrative signals not resolved by parcellation.")
  } else {
    report_lines <- c(report_lines, "",
                      "  Overall: Channel and source provide equivalent GAMM fit.",
                      "  Both domains appear to capture the alpha-pain relationship.")
  }
}

writeLines(report_lines, file.path(compare_out_dir, "comparison_report.txt"))

message("\nComparison complete. Outputs saved to: ", compare_out_dir)
message(readLines(file.path(compare_out_dir, "comparison_report.txt")), sep = "\n")