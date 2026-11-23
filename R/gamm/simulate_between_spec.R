# =============================================================================
# simulate_between_spec.R
# -----------------------------------------------------------------------------
# Which specification of the BETWEEN-subject (trait) term can detect a real
# between-person effect, while keeping false positives at ~5% when there is
# none? Decided on simulated data with the structure of the real screen, BEFORE
# the real data are refit (so the choice cannot be steered by real p-values).
#
# Why: in the first S1-cluster screen, 123 of 132 between-subject smooths had
# edf ~ 0 and no p < 0.05 (6.6 expected by chance). With ~57 trials per person
# the subject random intercept is almost unpenalized (edf 84 of 86) and
# absorbs the between-person variance; with select = TRUE the metric's
# subject-mean smooth is penalized INCLUDING its linear part, so fREML gives
# the variance to the random intercept and shrinks the metric term to zero.
#
# Specifications compared (same base model as gamm_screen.R):
#   A  current   s(x_between) + s(x_within), select = TRUE
#   B            s(x_between) + s(x_within), select = FALSE
#                (linear part unpenalized; wiggliness still penalized)
#   C            x_between (linear fixed effect, t-test) + s(x_within),
#                select = TRUE
# Truths: none (false-positive rate), linear r = 0.20 and 0.30 between the
# trait and the person's pain level, and a U-shaped effect (r = 0.30).
#
# Data shaped like the screen: 5 experiments (3 cap sizes), 18 subjects each
# (90), ~57 trials per subject, subject SD ~ residual SD (as in the fitted
# models: scale ~3.2, subject RE edf 84/86), laser power strong, metric =
# subject trait + experiment shift + trial noise (trial-level ICC ~0.5).
#
# Run from R/gamm/:  Rscript simulate_between_spec.R [n_sim] [n_cores] [truths] [specs] [seed_base]
#   defaults: n_sim = 200 per condition, cores = all but one (fork; Linux/mac),
#   all truths and specs, seed_base 0. Seeds are fixed so a run is reproducible:
#   the same arguments give the same datasets. Use another seed_base (e.g.
#   100000) for an INDEPENDENT batch.
#
# Results (Nov 2026; 90 subjects x 57 trials, 5 experiments):
#   False positives (no effect), 1000 datasets: A 8.1%, C 6.9% (nominal 5%;
#   the seed_base = 0 batch of 200 gave 11% for A-C, the 800-dataset
#   independent batch 7.4% / 6.0%). For reference on the same data, lmer
#   4.5% and a subject-level regression 3.5% (200 datasets).
#   Power at p < 0.05 (200 each): r = 0.20: A 46%, B 46%, C 45%;
#   r = 0.30: A 77%, B 78%, C 78%; U-shape r = 0.30: A 40%, B 37%, C 15%.
#   Conclusion: keep A (the screen as is). Its between-subject p-values are
#   mildly liberal (~7-8% at a nominal 5%) - acceptable for an uncorrected
#   screen whose hits must replicate on held-out subjects; report it.
#   Runtime ~ (4 truths x n_sim x 3 specs) fits of ~0.3-1 s each.
# Writes: simulate_between_spec_results.csv (one row per sim x spec) and
#         prints the rejection rate per truth x spec.
# =============================================================================

suppressPackageStartupMessages({ library(mgcv); library(parallel) })

args    <- commandArgs(trailingOnly = TRUE)
n_sim   <- if (length(args) >= 1) as.integer(args[1]) else 200L
n_cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, detectCores() - 1L)
sel_truths <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else NULL
sel_specs  <- if (length(args) >= 4) strsplit(args[4], ",")[[1]] else NULL
seed_base  <- if (length(args) >= 5) as.integer(args[5]) else 0L   # change for an independent batch
alpha   <- 0.05

sim_data <- function(effect, r = 0, n_exp = 5, n_per = 18, n_tr = 57) {
  n_s <- n_exp * n_per
  exp_id <- rep(seq_len(n_exp), each = n_per)
  cap <- c("62", "32", "32", "64", "32")[exp_id]
  trait <- rnorm(n_s)
  sd_subj <- 1.8                                      # total between-person SD of pain level
  shape <- switch(effect, none = 0 * trait, linear = trait,
                  ushape = (trait^2 - 1) / sqrt(2))   # unit variance either way
  b <- if (effect == "none") 0 else r * sd_subj
  u <- rnorm(n_s, 0, sqrt(sd_subj^2 - b^2))           # rest of the person's pain level
  e_shift_y <- rnorm(n_exp, 0, 0.3)
  e_shift_x <- rnorm(n_exp, 0, 0.5)                   # hardware / lab differences in the metric
  s <- data.frame(g = factor(seq_len(n_s)), e = factor(exp_id), cap = factor(cap),
                  age = rnorm(n_s), sex = factor(sample(c("F", "M"), n_s, TRUE)),
                  trait = trait, mu = 6 + b * shape + u + e_shift_y[exp_id],
                  xs = trait + e_shift_x[exp_id])
  d <- s[rep(seq_len(n_s), each = n_tr), ]
  d$lp <- rnorm(nrow(d))
  d$ti <- rep(seq_len(n_tr), n_s)
  d$x  <- d$xs + rnorm(nrow(d), 0, 1)                 # trial-level ICC of the metric ~0.5
  d$y  <- d$mu + 1.2 * d$lp + rnorm(nrow(d), 0, 1.8)
  d$xz <- as.numeric(scale(d$x))
  d$xb <- ave(d$xz, d$g)
  d$xw <- d$xz - d$xb
  d$lpz <- as.numeric(scale(d$lp)); d$tiz <- as.numeric(scale(d$ti))
  d
}

base <- "y ~ s(lpz, k = 8) + s(tiz, k = 10) + age + sex + cap + s(g, bs = 're') + s(e, bs = 're')"
specs <- list(
  A_current_select_true = list(f = paste(base, "+ s(xb, k = 10) + s(xw, k = 10)"), select = TRUE,  term = "s(xb)"),
  B_select_false        = list(f = paste(base, "+ s(xb, k = 10) + s(xw, k = 10)"), select = FALSE, term = "s(xb)"),
  C_linear_between      = list(f = paste(base, "+ xb + s(xw, k = 10)"),            select = TRUE,  term = "xb")
)
truths <- list(none = list(effect = "none", r = 0), linear_r0.20 = list(effect = "linear", r = 0.20),
               linear_r0.30 = list(effect = "linear", r = 0.30), ushape_r0.30 = list(effect = "ushape", r = 0.30))

if (!is.null(sel_truths)) truths <- truths[sel_truths]
if (!is.null(sel_specs))  specs  <- specs[sel_specs]

one_rep <- function(i, truth) {
  set.seed(seed_base + 1000 * match(truth$effect, c("none", "linear", "ushape")) + round(100 * truth$r) * 7919 + i)
  d <- sim_data(truth$effect, truth$r)
  do.call(rbind, lapply(names(specs), function(sn) {
    sp <- specs[[sn]]
    m <- tryCatch(bam(as.formula(sp$f), data = d, method = "fREML", select = sp$select,
                      discrete = FALSE, nthreads = 1), error = function(e) NULL)
    if (is.null(m)) return(data.frame(spec = sn, p = NA, edf = NA))
    sm <- summary(m)
    if (sp$term %in% rownames(sm$s.table)) {
      data.frame(spec = sn, p = sm$s.table[sp$term, "p-value"], edf = sm$s.table[sp$term, "edf"])
    } else {
      data.frame(spec = sn, p = sm$p.table[sp$term, 4], edf = 1)
    }
  }))
}

t0 <- Sys.time()
res <- do.call(rbind, lapply(names(truths), function(tn) {
  cat(sprintf("Truth %-14s %d sims x %d specs ...\n", tn, n_sim, length(specs)))
  out <- mclapply(seq_len(n_sim), one_rep, truth = truths[[tn]], mc.cores = n_cores)
  r <- do.call(rbind, Map(function(x, i) cbind(sim = i, x), out, seq_along(out)))
  cbind(truth = tn, r)
}))
write.csv(res, "simulate_between_spec_results.csv", row.names = FALSE)

cat(sprintf("\nRejection rate at p < %.2f (n_sim = %d per truth; MC error ~ +/-%.1f pts at 5%%, +/-%.1f at 50%%)\n",
            alpha, n_sim, 100 * 1.96 * sqrt(.05 * .95 / n_sim), 100 * 1.96 * sqrt(.25 / n_sim)))
tab <- aggregate(p ~ truth + spec, res, function(p) round(100 * mean(p < alpha), 1))
wide <- reshape(tab, idvar = "spec", timevar = "truth", direction = "wide")
names(wide) <- sub("^p\\.", "", names(wide))
print(wide[, c("spec", intersect(names(truths), names(wide)))], row.names = FALSE)
cat("\nShare of fits with the between term shrunk to zero (edf < 0.01), %:\n")
print(reshape(aggregate(edf ~ truth + spec, res, function(e) round(100 * mean(e < 0.01), 1)),
              idvar = "spec", timevar = "truth", direction = "wide"), row.names = FALSE)
cat("\nShare of p > 0.90, % (real screen: 23 of 132 = 17%):\n")
print(reshape(aggregate(p ~ truth + spec, res, function(p) round(100 * mean(p > 0.9), 1)),
              idvar = "spec", timevar = "truth", direction = "wide"), row.names = FALSE)
cat("\nMedian edf of the between term:\n")
print(reshape(aggregate(edf ~ truth + spec, res, function(e) round(median(e), 3)),
              idvar = "spec", timevar = "truth", direction = "wide"), row.names = FALSE)
cat(sprintf("\nfailed fits: %d | runtime %.1f min\n", sum(is.na(res$p)),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("\n'none' column = false-positive rate (should be ~5%). Other columns = power.\n")