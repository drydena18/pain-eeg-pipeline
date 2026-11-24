# =============================================================================
# simulate_within_holdout_power.R                                     V1.0.0
# -----------------------------------------------------------------------------
# Power of the two within-person held-out criteria (holdout_replicate.R V4)
# for a given within-person effect, by simulation. Run from R/gamm/:
#   Rscript simulate_within_holdout_power.R <beta> <n_sim> <seed> [n_per_exp] [n_holdout_per_exp]
# beta = SD of the true within-person effect in rating points (0.065 ~ a 0.21-
# point 5th-95th percentile spread, the size seen for S1_contra pre_sf_logratio).
# 4 experiments, 53 trials/subject, residual SD 1.78 (matching the Nov 2026 screen).
#
# Results (Nov 2026, 0-10 scale):
#   beta 0.065, 40 held-out : screen 78 %, prediction 25 %, refit 38 %
#   beta 0.10,  40 held-out : screen 100 %, prediction 28 %, refit 88 %
#   beta 0,     40 held-out : screen 12 %, prediction 5 %, refit 7 %  (false positives)
#   beta 0.065, 128 held-out: prediction 40-47 %, refit 87-100 %
# =============================================================================
suppressPackageStartupMessages({library(mgcv); library(dplyr)})
source("gamm_defaults.R"); source("helpers/gamm_combo_core.R"); source("helpers/hit_helpers.R")
args <- commandArgs(TRUE); beta <- as.numeric(args[1]); nsim <- as.integer(args[2]); seed0 <- as.integer(args[3])
n_exp_subj <- if (length(args) >= 4) as.integer(args[4]) else 40L
n_ho_subj  <- if (length(args) >= 5) as.integer(args[5]) else 10L
cfg <- gamm_cfg; cfg$nthreads <- 1
one <- function(seed) {
  set.seed(seed); nS <- 4 * n_exp_subj; nT <- 53
  subj <- tibble(global_subjid = sprintf("E%02d_S%03d", rep(1:4, each = n_exp_subj), rep(seq_len(n_exp_subj), 4)),
                 experiment_id = rep(as.character(1:4), each = n_exp_subj),
                 split = rep(c(rep("train", n_exp_subj - n_ho_subj), rep("holdout", n_ho_subj)), 4),
                 age = round(rnorm(nS, 25, 4)), sex = sample(c("F", "M"), nS, TRUE),
                 cap_size = rep(c("32", "32", "64", "64"), each = n_exp_subj), u = rnorm(nS, 0, 1.2), xm = rnorm(nS, 0, 0.6))
  d <- subj[rep(seq_len(nS), each = nT), ] %>% group_by(global_subjid) %>%
    mutate(trial_index = row_number(), laser_power = sample(rep(seq(2, 4.5, by = 0.25), length.out = nT))) %>% ungroup() %>%
    mutate(x = xm + rnorm(n(), 0, 0.5), pain_rating = 4 + 1.6 * laser_power + u + beta * (x - xm) / 0.5 + rnorm(n(), 0, 1.78))
  zs <- list(); zf <- function(x, v, s = NULL) { if (is.null(s)) { zs[[v]] <<- c(mean = mean(x), sd = sd(x)); s <- zs[[v]] }; (x - s[["mean"]]) / s[["sd"]] }
  prep <- function(df, train) df %>% mutate(age_z = zf(age, "age", if (!train) zs$age),
      laser_power_z = zf(laser_power, "laser_power", if (!train) zs$laser_power),
      trial_index_z = zf(trial_index, "trial_index", if (!train) zs$trial_index),
      x_z = zf(x, "x", if (!train) zs$x)) %>% add_between_within("x_z")
  tr <- prep(filter(d, split == "train"), TRUE) %>% mutate(across(c(global_subjid, experiment_id, sex, cap_size), factor))
  ho <- prep(filter(d, split == "holdout"), FALSE)
  f <- pain_rating ~ s(laser_power_z, k = 10) + s(trial_index_z, k = 10) + age_z + sex + cap_size +
    s(global_subjid, bs = "re") + s(experiment_id, bs = "re") + s(x_between_z, k = 10) + s(x_within_z, k = 10)
  m <- bam(f, data = tr, method = "fREML", select = TRUE)
  m$gamm_meta <- list(decomposition = "between_within", primary_part = "between")
  screen_p <- summary(m)$s.table["s(x_within_z)", "p-value"]
  base <- refit_baseline(m, "x_z", cfg, drop = "x_within_z")
  ho$pred_metric <- predict_population(m, ho); ho$pred_base <- predict_population(base, ho)
  we <- within_person_errors(ho); bw <- boot_subjects(we$num, we$den, 500)
  rf <- refit_holdout(m, ho, "x_within_z", cfg)
  c(screen = screen_p < .05, pred = unname(bw[["p"]] < .05 && bw[["est"]] < 0),
    refit = rf$p < .05 && rf$direction > 0, pct = 100 * bw[["est"]] / (sum(we$sse_b) / sum(we$den)))
}
r <- t(sapply(seed0 + seq_len(nsim), one))
cat(sprintf("beta = %.3f (5-95%% spread ~ %.2f points): screen power %.2f | held-out prediction %.2f | held-out refit %.2f | mean pct %.3f\n",
            beta, beta * 3.29, mean(r[, "screen"]), mean(r[, "pred"]), mean(r[, "refit"]), mean(r[, "pct"])))