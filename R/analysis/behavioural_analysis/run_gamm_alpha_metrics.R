# ============================================================
# run_gamm_alpha_metrics.R
# ============================================================
# Consolidates trial-level GAMM pipeline for alpha-pain analyses
# across all experiments.
#
# Model set (all models include the same random-effects and
# covariate base):
#   m00 baseline (covariates only)
#   m01 baseline + aperiodic (fooof offset + exponent)
#   m02 sf_logratio
#   m03 sf_balance
#   m04 slow_frac
#   m05 te(pow_slow, pow_fast) [tensor interaction]
#   m06 paf_cog_hz
#   m07 sf_logratio + aperiodic
#   m08 sf_balance + aperiodic
#   m09 slow_frac + aperiodic
#   m10 te(pow_slow, pow_fast) + aperiodic
#   m11 pow_slow + pow_fast (main effects only) + aperiodic
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(mgcv)
library(tibble)
library(ggplot2)

# ============================================================
# USER SETTINGS
# ============================================================
data_file         <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/alpha_metrics_data.csv"
out_dir           <- "pain-eeg-pipeline/R/analysis/behavioural_analysis/gamm_outputs"
min_fooof_r2      <- 0.80
max_fooof_error   <- Inf 
use_discrete      <- TRUE
nthreads_to_use   <- 4L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# HELPERS
# ============================================================
clean_names_local <- function(x) {
    x %>%
        str_trim() %>%
        str_replace_all("\\s+", "_") %>%
        str_replace_all("\\^2", "r2") %>%
        str_replace_all("\\.", "_")
}

zscore <- function(x) as.numeric(scale(x))
safe_factor <- function(x) as.factor(as.character(x))

# Add a single smooth term to an existing formula
add_smooth <- function(base_f, var_z, k = 10) {
    new_term <- paste0("s(", var_z, ", k = ", k, ")")
    old_rhs <- deparse(base_f[[3L]])
    new_rhs <- paste(old_rhs, "+", new_term)
    as.formula(paste(deparse(base_f[[2L]]), "~", new_rhs))
}

# Add a tensor smooth to an existing formula
add_tensor <- function(base_f, var1, var2, k = 10) {
    new_term <- paste0("te(", var1, ", ", var2, ", k = ", k, ")")
    old_rhs <- deparse(base_f[[3L]])
    new_rhs <- paste(old_rhs, "+", new_term)
    as.formula(paste(deparse(base_f[[2L]]), "~", new_rhs))
}

extract_model_info <- function(mod, nm) {
    tibble(
        model_name = nm,
        AIC = AIC(mod),
        BIC = BIC(mod),
        logLik = as.numeric(logLik(mod)),
        dev_expl = summary(mod)$dev.expl,
        n = nobs(mod)
    )
}

plot_smooth_term <- function(model, term_z, file_name, out_dir) {
    labels <- sapply(model$smooth, function(s) s$label)
    pattern <- paste0("(^|,\\s*)", term_z, "(\\s*,|\\))")
    sel <- which(str_detect(labels, pattern))

    if (length(sel) == 0L) {
        warning("Term '", term_z, "' not found in model smooths. Skipping plot.")
        return(invisible(NULL))
    }

    png(file.path(out_dir, file_name), width = 1600, height = 1200, res = 200)
    plot(
        model,
        select = sel[1L],
        shade = TRUE,
        shade.col = "lightblue",
        main = paste("Effect of", term_z, "on Pain Rating."),
        xlab = term_z,
        ylab = "Partial Effect on Pain"
    )
    abline(h = 0, lty = 2)
    dev.off()
}

fit_bam <- function(formula, data, label) {
    message("Fitting: ", label)
    tryCatch(
        bam(formula = formula,
        data = data,
        method = 'fREML',
        discrete = use_discrete,
        nthreads = nthreads_to_use),
        error = function(e) {
            warning("Model '", label, "' failed: ", conditionMessage(e))
            NULL
        }
    )
}

# ============================================================
# LOAD DATA
# ============================================================
if (!file.exists(data_file)) stop("Data file not found: ", data_file)

df <- read_csv(data_file, show_col_types = FALSE)
names(df) <- clean_names_local(names(df))

required_cols <- c(
    "experiment_name", "experiment_id", "subjid", "subjid_uid", "global_subjid",
    "age", "sex", "cap_size", "trial", "trial_index", "laser_power", "pain_rating"
)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
    mutate(
        experiment_name = safe_factor(experiment_name),
        experiment_id = safe_factor(experiment_id),
        subjid = as.integer(subjid),
        subjid_uid = safe_factor(subjid_uid),
        global_subjid = safe_factor(global_subjid),
        age = suppressWarnings(as.numeric(age)),
        sex = safe_factor(sex),
        cap_size = safe_factor(cap_size),
        trial = as.integer(trial),
        trial_index = as.integer(trial_index),
        laser_power = suppressWarnings(as.numeric(laser_power)),
        pain_rating = suppressWarnings(as.numeric(pain_rating))
    )

# ============================================================
# QC FILTERING
# ============================================================
alpha_metrics <- c("pow_slow_alpha", "pow_fast_alpha", "sf_logratio", "sf_balance", "slow_alpha_frac")
present_alpha <- intersect(alpha_metrics, names(df))

df_model <- df %>%
    filter(
        !is.na(pain_rating),
        !is.na(laser_power),
        !is.na(trial_index),
        !is.na(global_subjid),
        !is.na(experiment_id),
        !is.na(age),
        !is.na(sex),
        !is.na(cap_size)
    )

# Filter out rows missing ALL alpha metrics
if (length(present_alpha) > 0) {
    df_model <- df_model %>%
        filter(if_any(all_of(present_alpha), ~ !is.na(.)))
}

if ("fooof_r2" %in% names(df_model)) {
    df_model <- df_model %>%
        filter(is.na(fooof_r2) | fooof_r2 >= min_fooof_r2)
}

if ("fooof_error" %in% names(df_model) && is.finite(max_fooof_error)) {
    df_model <- df_model %>%
        filter(is.na(fooof_error) | fooof_error <= max_fooof_error)
}

message("Rows after QC filtering: ", nrow(df_model))

# ============================================================
# Z-SCORE CONTINUOUS PREDICTORS
# ============================================================
# Core predictors (always present)
df_model <- df_model %>%
    mutate(
        age_z = zscore(age),
        laser_power_z = zscore(laser_power),
        trial_index_z = zscore(trial_index)
    )

# Optional alpha metrics - z-score only if present
optional_metrics <- c(
    "paf_cog_hz", "pow_slow_alpha", "pow_fast_alpha",
    "sf_ratio", "sf_balance", "sf_logratio",
    "slow_alpha_frac", "rel_slow_alpha", "rel_fast_alpha",
    "fooof_offset", "fooof_exponent"
)
for (m in optional_metrics) {
    if (m %in% names(df_model)) {
        z_col <- paste0(m, "_z")
        df_model[[z_col]] <- zscore(df_model[[m]])
    }
}

write_csv(df_model, file.path(out_dir, "alpha_pain_master_model_input.csv"))
message("Model input saved: ", file.path(out_dir, "alpha_pain_master_model_input.csv"))

# ============================================================
# CHECK COLUMN AVAILABILITY
# ============================================================
has_col <- function(col) col %in% names(df_model)

has_slow <- has_col("pow_slow_alpha_z")
has_fast <- has_col("pow_fast_alpha_z")
has_sf_logratio <- has_col("sf_logratio_z")
has_sf_balance <- has_col("sf_balance_z")
has_slow_frac <- has_col("slow_alpha_frac_z")
has_paf <- has_col("paf_cog_hz_z")
has_fooof <- has_col("fooof_offset_z") && has_col("fooof_exponent_z")

# ============================================================
# BASE FORMULAS
# ============================================================
base_formula <- pain_rating ~
    s(laser_power_z, k = 10) +
    s(trial_index_z, k = 10) +
    s(age_z, k = 10) +
    sex +
    cap_size +
    s(global_subjid, bs = "re") +
    s(experiment_id, bs = "re")

# Build aperiodic base formula
if (has_aperiodic) {
    base_formula_aperiodic <- pain_rating ~
        s(laser_power_z, k = 10) +
        s(trial_index_z, k = 10) +
        s(age_z, k = 10) +
        sex +
        cap_size +
        s(fooof_offset_z, k = 10) +
        s(fooof_exponent_z, k = 10) +
        s(global_subjid, bs = "re") +
        s(experiment_id, bs = "re")
} else {
    message("[WARN] FOOOF columns (fooof_offset_z, fooof_exponent_z) not found. "
            "Aperiodic models (m01, m07-m11) will be skipped.")
    base_formula_aperiodic <- NULL
}

# ============================================================
# FIT MODELS
# ============================================================
model_list <- list()

model_list[["m00_baseline"]] <- fit_bam(base_formula, df_model, "m00_baseline")

if (!is.null(base_formula_aperiodic)) {
    model_list[["m01_baseline_aperiodic"]] <-
        fit_bam(base_formula_aperiodic, df_model, "m01_baseline_aperiodic")
}

if (has_sflog) {
    model_list[["m02_sf_logratio"]] <-
        fit_bam(add_smooth(base_formula, "sf_logratio_z"), df_model, "m02_sf_logratio")
}

if (has_sfbal) {
    model_list[["m03_sf_balance"]] <-
        fit_bam(add_smooth(base_formula, "sf_balance_z"), df_model, "m03_sf_balance")
}

if (has_slowfrac) {
    model_list[["m04_slow_frac"]] <-
        fit_bam(add_smooth(base_formula, "slow_alpha_frac_z"), df_model, "m04_slow_frac")
}

if (has_slow && has_fast) {
    model_list[["m05_pow_slow_fast_te"]] <-
        fit_bam(add_tensor(base_formula, "pow_slow_alpha_z", "pow_fast_alpha_z"), df_model, "m05_pow_slow_fast_te")
}

if (has_paf) {
    model_list[["m06_paf_cog_hz"]] <-
        fit_bam(add_smooth(base_formula, "paf_cog_hz_z"), df_model, "m06_paf_cog_hz")
}

if (!is.null(base_formula_aperiodic)) {
    if (has_sflog) {
        model_list[["m07_sf_logratio_aperiodic"]] <-
            fit_bam(add_smooth(base_formula_aperiodic, "sf_logratio_z"), df_model, "m07_sf_logratio_aperiodic")
    }
    if (has_sfbal) {
        model_list[["m08_sf_balance_aperiodic"]] <-
            fit_bam(add_smooth(base_formula_aperiodic, "sf_balance_z"), df_model, "m08_sf_balance_aperiodic")
    }
    if (has_slowfrac) {
        model_list[["m09_slow_frac_aperiodic"]] <-
            fit_bam(add_smooth(base_formula_aperiodic, "slow_alpha_frac_z"), df_model, "m09_slow_frac_aperiodic")
    }
    if (has_slow && has_fast) {
        model_list[["m10_pow_slow_fast_te_aperiodic"]] <-
            fit_bam(add_tensor(base_formula_aperiodic, "pow_slow_alpha_z", "pow_fast_alpha_z"), df_model, "m10_pow_slow_fast_te_aperiodic")
        model_list[["m11_pow_slow_fast_main_aperiodic"]] <-
            fit_bam(
                update(base_formula_aperiodic, . ~ . + pow_slow_alpha_z + pow_fast_alpha_z),
                df_model,
                "m11_pow_slow_fast_main_aperiodic"
            )
    }
}

# Drop any models that failed to fit
model_list <- compact(model_list)

if (length(model_list) == 0L) stop("No GAMMs were successfully fit.")
message(length(model_list), " models fit successfully.")

# ============================================================
# SAVE SUMMARIES
# ============================================================
iwalk(model_list, function(mod, nm) {
    capture.output(
        summary(mod),
        file = file.path(out_dir, paste0(nm, "_summary.txt"))
    )
})

# ============================================================
# MODEL COMPARISON TABLE
# ============================================================
model_comparison <- imap_dfr(model_list, extract_model_info) %>%
    arrange(AIC) %>%
    mutate(delta_AIC = AIC - min(AIC))

write_csv(model_comparison, file.path(out_dir, "model_comparison.csv"))
message("Model comparison saved.")
print(model_comparison)

# ============================================================
# SMOOTH EFFECTS PLOTS
# ============================================================
smooth_plots <- list(
    m02_sf_logratio = "sf_logratio_z",
    m03_sf_balance = "sf_balance_z",
    m04_slow_frac = "slow_alpha_frac_z",
    m06_paf_cog_hz = "paf_cog_hz_z",
    m07_sf_logratio_aperiodic = "sf_logratio_z",
    m08_sf_balance_aperiodic = "sf_balance_z",
    m09_slow_frac_aperiodic = "slow_alpha_frac_z"
)

for (nm in names(smooth_plots)) {
    if (nm %in% names(model_list)) {
        plot_smooth_term(
            model = model_list[[nm]],
            term_z = smooth_plots[[nm]],
            file_name = paste0(nm, "_effect.png"),
            out_dir = out_dir
        )
    }
}

# ============================================================
# INTERACTION PLOTS (if tensor terms were fit)
# ============================================================
tensor_models <- c("m05_pow_slow_fast_te", "m10_pow_slow_fast_te_aperiodic")

for (nm in tensor_models) {
    if (nm %in% names(model_list)) {
        png(file.path(out_dir, paste0(nm, "_surface.png")),
        width = 1800, height = 1400, res = 200)
    vis.gam(
        model_list[[nm]],
        view = c("pow_slow_alpha_z", "pow_fast_alpha_z"),
        plot.type = "contour",
        too.far = 0.05,
        main = paste("Slow-Fast Alpha Tensor:", nm)
    )
    dev.off()
    }
}

# ============================================================
# FITTED VALUES / RESIDUALS
# ============================================================
# Selected output columns that actually exist
fitted_base_cols <- c(
    "experiment_name", "experiment_id",
    "subjid", "subjid_uid", "global_subjid",
    "age", "sex", "cap_size",
    "trial", "trial_index", "laser_power", "pain_rating"
)
fitted_optional_cols <- c(
    "paf_cog_hz", "pow_slow_alpha", "pow_fast_alpha",
    "sf_ratio", "sf_balance", "sf_logratio",
    "slow_alpha_frac", "fooof_offset", "fooof_exponent"
)

fitted_trials_df <- df_model %>%
    select(all_of(fitted_base_cols), any_of(fitted_optional_cols))

for (nm in names(model_list)) {
    fitted_trial_df[[paste0(nm, "_fitted")]] <- fitted(model_list[[nm]])
    fitted_trial_df[[paste0(nm, "_resid")]] <- residuals(model_list[[nm]])
}

write_csv(fitted_trial_df, file.path(out_dir, "trial_level_fitted_values.csv"))
message("Trial-level fitted values saved.")

# Subject-level grand average summary
subject_ga_df <- fitted_trial_df %>%
    group_by(
        experiment_name, experiment_id,
        subjid, subjid_uid, global_subjid,
        age, sex, cap_size
    ) %>%
    summarise(
        n_trials = n(),
        mean_pain_rating = mean(pain_rating, na.rm = TRUE),
        mean_laser_power = mean(laser_power, na.rm = TRUE),
        across(any_of(fitted_optional_cols), ~ mean(.x, na.rm = TRUE), .names = "mean_{col}"),
        across(ends_with("_fitted"), ~ mean(.x, na.rm = TRUE), .names = "ga_{col}"),
        .groups = "drop"
    )

write_csv(subject_ga_df, file.path(out_dir, "subject_level_ga_fitted_values.csv"))
message("Subject-level GA fitted values saved.")

# ============================================================
# DIAGNOSTICS
# ============================================================
iwalk(model_list, function(mod, nm) {
    png(file.path(out_dir, paste0(nm, "_diagnostics.png")),
    width = 1800, height = 1400, res = 180)
    gam.check(mod)
    dev.off()
})

# Observed v Fitted scatter plts
for (nm in names(model_list)) {
    p <- ggplot(
        fitted_trial_df,
        aes(x = .data[[paste0(nm, "_fitted")]], y = pain_rating)
    ) +
    geom_point(alpha = 0.15, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
    labs(
        title = paste("Observed vs Fitted:", nm),
        x = "Fitted pain rating",
        y = "Observed pain rating"
    ) +
    theme_minimal()

ggsave(
    filename = file.path(out_dir, paste0(nm, "_observed_vs_fitted.png")),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
)
}

# ============================================================
# SAVE MODEL OBJECTS
# ============================================================
iwalk(model_list, function(mod, nm) {
    saveRDS(mod, file.path(out_dir, paste0(nm, ".rds")))
})

message("GAMM pipeline complete. Output: ", out_dir)