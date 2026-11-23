# =============================================================================
# test_missing_covariates.R
# -----------------------------------------------------------------------------
# Checks gamm_cfg$missing_covariates = "impute" (Nov 2026):
#   impute_covariates() / apply_covariate_fill() / adds_rank()
#   (helpers/gamm_combo_core.R), a real fit_gamm_combo() fit on simulated
#   subjects with missing age / sex, and held-out prediction for people with
#   missing age / sex filled with the TRAINING values (as prepare_holdout()
#   does). No database needed.
#
# Run from R/gamm/:  Rscript ../tests/test_missing_covariates.R
# Needs mgcv, dplyr, tibble, purrr, stringr.
# =============================================================================
suppressPackageStartupMessages({ library(mgcv); library(dplyr); library(tibble); library(purrr); library(stringr) })

load_funs <- function(file, names) for (e in parse(file))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% names) eval(e, globalenv())
load_funs("helpers/gamm_combo_core.R",
          c("%||%", "MISSING_COVARIATES", "impute_covariates", "apply_covariate_fill", "adds_rank",
            "safe_k", "smooth_or_linear", "safe_bam", "metric_parts", "add_between_within",
            "build_combo_formula", "smooth_term_variables", "extract_combo_result_row", "fit_gamm_combo"))
load_funs("helpers/hit_helpers.R", "predict_population")

n_pass <- 0; n_all <- 0
check <- function(cond, msg) {
  n_all <<- n_all + 1
  if (isTRUE(cond)) { n_pass <<- n_pass + 1; cat("PASS ", msg, "\n") } else cat("FAIL ", msg, "\n")
}

# -- Simulated training data: 40 subjects, 3 experiments; some missing age / sex ----
set.seed(7)
mk <- function(n, ids, tr = 30) {
  s <- data.frame(global_subjid = ids, experiment_id = rep(c("1", "2", "3"), length.out = n),
                  age = round(rnorm(n, 25, 4)), sex = sample(c("F", "M"), n, TRUE))
  s$cap_size <- c("1" = "62", "2" = "32", "3" = "64")[s$experiment_id]
  d <- s[rep(seq_len(n), times = c(rep(tr, n - 1), 3 * tr)), ]    # last subject has 3x the trials
  d$trial_index <- ave(seq_len(nrow(d)), d$global_subjid, FUN = seq_along)
  d$laser_power <- rnorm(nrow(d)); d$m <- rnorm(nrow(d))
  d$pain_rating <- 5 + d$laser_power + rep(rnorm(n), table(factor(d$global_subjid, ids))) + rnorm(nrow(d))
  d
}
ids <- sprintf("S%02d", 1:40)
tr <- mk(40, ids)
tr$age[tr$global_subjid %in% c("S01", "S02", "S03")] <- NA         # 3 missing age
tr$sex[tr$global_subjid %in% c("S04", "S05")] <- NA                # 2 missing sex
tr$sex[tr$global_subjid == "S06"] <- ""                            # 1 empty-string sex

# -- impute_covariates ----
dr <- impute_covariates(tr, "drop")
check(is.null(dr$fill) && identical(dr$df, tr), "drop: data unchanged, no fill")
im <- impute_covariates(tr, "impute")
subj_age <- tr$age[!duplicated(tr$global_subjid)]
check(im$fill$n_age == 3 && im$fill$n_sex == 3, "impute: counts subjects (3 age, 3 sex incl. empty string)")
check(im$fill$age == median(subj_age, na.rm = TRUE), "impute: age = median over SUBJECTS (not trials)")
check(!any(is.na(im$df$age)) && all(im$df$sex %in% c("F", "M", "U")), "impute: no NA age / sex left; unknown sex = 'U'")
check(sum(im$df$age_missing[!duplicated(im$df$global_subjid)]) == 3, "impute: age_missing = 1 for exactly the 3 subjects")
err <- tryCatch({ impute_covariates(tr, "mean"); FALSE }, error = function(e) TRUE)
check(err, "unknown missing_covariates value -> error")

# -- adds_rank ----
d1 <- im$df; d1$cap_size <- factor(d1$cap_size); d1$sex <- factor(d1$sex)
check(adds_rank(c("sex", "cap_size"), "age_missing", d1), "age_missing added when only some subjects lack age")
d2 <- d1; d2$age_missing <- as.numeric(d2$cap_size == "64")          # age missing for one whole cap size
check(!adds_rank(c("sex", "cap_size"), "age_missing", d2), "age_missing skipped when it duplicates a cap size")

# -- Fit as the screen does (gamm_screen.R logic) ----
df <- im$df %>% mutate(global_subjid = factor(global_subjid), experiment_id = factor(experiment_id),
                       sex = factor(sex), cap_size = factor(cap_size),
                       age_z = as.numeric(scale(age)), laser_power_z = as.numeric(scale(laser_power)),
                       trial_index_z = as.numeric(scale(trial_index)), m_z = as.numeric(scale(m)))
zscale <- lapply(c(age = "age", laser_power = "laser_power", trial_index = "trial_index", m = "m"),
                 function(cn) c(mean = mean(df[[cn]]), sd = sd(df[[cn]])))
df <- add_between_within(df, "m_z")
fixed <- c("age_z", "sex", "cap_size")
if (adds_rank(fixed, "age_missing", df)) fixed <- c(fixed, "age_missing")
base <- as.formula(paste("pain_rating ~ s(laser_power_z, k = 8) + s(trial_index_z, k = 10) +",
                         paste(c(fixed, "s(global_subjid, bs = 're')", "s(experiment_id, bs = 're')"), collapse = " + ")))
out_dir <- tempfile("gamm_test_"); dir.create(out_dir)
fit <- fit_gamm_combo(df, "channel", "S1_contra", "m_z", base, "m_1metric_m", out_dir,
                      zscale = zscale, decomposition = "between_within", nthreads = 1L,
                      covariate_fill = im$fill)
mod <- fit$model
check(!is.null(mod), "screen-style fit with imputed covariates succeeds")
check(n_distinct(mod$model$global_subjid) == 40, "all 40 subjects in the model (complete cases would keep 34)")
check(all(c("sexU", "age_missing") %in% rownames(summary(mod)$p.table)), "model has sexU and age_missing terms")
check(identical(readRDS(file.path(out_dir, "m_1metric_m.rds"))$covariate_fill, im$fill), "fill saved on the .rds")
check(identical(mod$gamm_meta$missing_covariates, "impute"), "gamm_meta records missing_covariates = impute")

# -- Held-out people with missing age / sex (prepare_holdout path) ----
ho <- mk(6, sprintf("H%02d", 1:6))
ho$age[ho$global_subjid == "H01"] <- NA
ho$sex[ho$global_subjid == "H02"] <- NA
ho <- apply_covariate_fill(ho, mod$covariate_fill)
check(all(ho$age[ho$global_subjid == "H01"] == im$fill$age), "held-out age filled with the TRAINING median")
check(all(ho$sex[ho$global_subjid == "H02"] == "U"), "held-out unknown sex -> 'U'")
z <- function(x, v) (x - zscale[[v]][["mean"]]) / zscale[[v]][["sd"]]
ho <- ho %>% mutate(age_z = z(age, "age"), laser_power_z = z(laser_power, "laser_power"),
                    trial_index_z = z(trial_index, "trial_index"), m_z = z(m, "m"))
ho <- add_between_within(ho, "m_z")
pr <- predict_population(mod, ho)
check(all(!is.na(pr)), "predictions for every held-out person, including missing age / sex")

# -- Complete cases ("drop") for comparison ----
dd <- tr %>% filter(!is.na(age), !is.na(sex), sex != "")
check(n_distinct(dd$global_subjid) == 34, "drop mode would keep 34 of 40 subjects")

unlink(out_dir, recursive = TRUE)
cat(sprintf("\n%d / %d checks passed\n", n_pass, n_all))
if (n_pass != n_all) quit(status = 1)