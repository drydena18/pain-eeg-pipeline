# =============================================================================
# check_discrete_between_p.R
# -----------------------------------------------------------------------------
# Does bam(discrete = TRUE) give valid p-values for a subject-level
# (between-subject) smooth on THIS mgcv version? Simulated data; compare the
# two p columns: with a real effect (beta = 0.5) both should be tiny; with
# beta = 0 both should be non-extreme. On mgcv 1.9-1 the discrete p-values
# were ~1 in every case, which is why gamm_defaults.R sets bam_discrete = FALSE.
# Runs in about 30 s. Rscript check_discrete_between_p.R
# =============================================================================
suppressMessages(library(mgcv)); set.seed(3)
sim <- function(n_exp, n_per, n_tr, beta) {
  s <- data.frame(g = factor(seq_len(n_exp * n_per)), e = factor(rep(seq_len(n_exp), each = n_per)),
                  trait = rnorm(n_exp * n_per), u = rnorm(n_exp * n_per), age = rnorm(n_exp * n_per),
                  sex = factor(sample(c("F","M"), n_exp*n_per, TRUE)))
  d <- s[rep(seq_len(nrow(s)), each = n_tr), ]
  d$x <- d$trait + rnorm(nrow(d), 0, 0.6); d$lp <- rnorm(nrow(d))
  d$y <- 1 + d$lp + beta * d$trait + d$u + rnorm(nrow(d))
  d$xz <- as.numeric(scale(d$x)); d$xb <- ave(d$xz, d$g); d$xw <- d$xz - d$xb
  d
}
for (cfg in list(c(3,22,25,0.5), c(9,15,30,0.5), c(9,40,30,0.5), c(9,40,30,0))) {
  d <- sim(cfg[1], cfg[2], cfg[3], cfg[4])
  f <- y ~ s(lp) + age + sex + s(g, bs="re") + s(e, bs="re") + s(xb) + s(xw)
  p <- sapply(c(TRUE, FALSE), function(disc) { m <- bam(f, data=d, method="fREML", select=TRUE, discrete=disc, nthreads=1)
       summary(m)$s.table["s(xb)", "p-value"] })
  cat(sprintf("exps=%d subj=%d trials/subj=%d beta=%.1f | between p: discrete=%.3g  non-discrete=%.3g\n", cfg[1], cfg[1]*cfg[2], cfg[3], cfg[4], p[1], p[2]))
}