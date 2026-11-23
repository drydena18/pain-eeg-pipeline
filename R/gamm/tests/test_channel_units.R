# =============================================================================
# test_channel_units.R
# -----------------------------------------------------------------------------
# Checks gamm_cfg$analysis_levels and channel_units / cluster_names (V4, Nov 2026):
#   apply_channel_units(), select_screen_units() (the screen),
#   filter_channel_units() and load_hits() (every follow-up script),
#   get_analysis_levels() / check_screen_level(), and the classical-test
#   level / unit selection.
# No real data needed: builds a mock gamm_fitted in an in-memory SQLite DB
# that holds an older all-electrode screen alongside a cluster screen.
#
# Run from R/gamm/:  Rscript ../tests/test_channel_units.R
# Needs DBI + RSQLite (test only), dplyr, tibble.
# =============================================================================
suppressPackageStartupMessages({ library(DBI); library(dplyr); library(tibble) })
`%||%` <- function(a, b) if (is.null(a)) b else a

# Load only the functions under test (db_helpers.R also loads duckdb, which
# this test does not need).
load_funs <- function(file, names) {
  ex <- parse(file)
  for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
                    as.character(e[[2]]) %in% names) eval(e, globalenv())
}
load_funs("helpers/db_helpers.R", c("CHANNEL_UNITS", "ANALYSIS_LEVELS", "apply_channel_units", "select_units",
                                    "get_analysis_levels"))
load_funs("helpers/gamm_screen.R", c("select_screen_units", "check_screen_level"))
load_funs("classical_tests.R", c("select_classical_units", "classical_level"))
load_funs("helpers/hit_helpers.R", c("load_hits", "filter_channel_units"))

n_pass <- 0; n_all <- 0
check <- function(cond, msg) {
  n_all <<- n_all + 1
  if (isTRUE(cond)) { n_pass <<- n_pass + 1; cat("PASS ", msg, "\n") } else cat("FAIL ", msg, "\n")
}

clusters <- c("S1_contra", "S1_ipsi", "SM_furman")
avail <- c("C3", "C4", "CP1", "Cz", "Pz", "S1_contra", "S1_ipsi", "SM_furman")

# -- apply_channel_units ----
check(identical(apply_channel_units(avail, "clusters", clusters, FALSE), clusters), "clusters -> cluster rows only")
check(identical(apply_channel_units(avail, "electrodes", clusters, FALSE), c("C3", "C4", "CP1", "Cz", "Pz")), "electrodes -> electrode rows only")
check(identical(apply_channel_units(avail, "all", clusters, FALSE), avail), "all -> everything")
check(identical(apply_channel_units(avail, "CLUSTERS", tolower(clusters), FALSE), clusters), "case-insensitive (setting and names)")
err <- tryCatch({ apply_channel_units(avail, "cluster", clusters); FALSE }, error = function(e) TRUE)
check(err, "typo in channel_units -> error")
msgs <- character(0)
invisible(withCallingHandlers(apply_channel_units(c(avail[-8], "Foo_bar"), "clusters", clusters, TRUE),
                    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }))
check(any(grepl("SM_furman", msgs)) && any(grepl("Foo_bar", msgs)), "warns: listed cluster absent / unlisted '_' label")

# -- the screen's unit selection ----
cfg <- list(channel_units = "clusters", cluster_names = clusters, channel_include = NULL, roi_include = NULL)
check(identical(select_screen_units(avail, "channel", cfg), clusters), "screen: clusters")
cfg$channel_include <- "S1_contra"
check(identical(select_screen_units(avail, "channel", cfg), "S1_contra"), "screen: channel_include narrows within clusters")
cfg$channel_include <- "^S1_"
check(identical(select_screen_units(avail, "channel", cfg), c("S1_contra", "S1_ipsi")), "screen: regex include within clusters")
cfg2 <- modifyList(cfg, list(channel_units = "electrodes", channel_include = NULL))
check(identical(select_screen_units(avail, "channel", cfg2), c("C3", "C4", "CP1", "Cz", "Pz")), "screen: electrodes")
check(identical(select_screen_units(c("S1", "ACC"), "source", cfg), c("S1", "ACC")), "screen: source level unaffected")
check(is.null(cfg$no_such) && identical(suppressMessages(select_screen_units(avail, "channel", list())), avail),
      "screen: missing setting = all (old configs keep working)")

# -- load_hits on a mock gamm_fitted ----
con <- dbConnect(RSQLite::SQLite(), ":memory:")
mk <- function(level, unit, p, model = "m_1metric_pre_sf_balance", n_metrics = 1L)
  tibble(feature_set = "v3.1", decomposition = "between_within", level = level, group_value = unit,
         model_name = model, metrics = "pre_sf_balance_z", n_metrics = n_metrics, primary_part = "between",
         edf = "2.1", f_stat = "5.0", p_value = as.character(p), edf_other = "1", f_other = "1",
         p_other = "0.5", aic = 1000, dev_expl = 0.3, r_sq = 0.3, n_obs = 500L, n_subjects = 40L,
         bam_discrete = 0L, rds_path = NA_character_)
fitted <- bind_rows(
  mk("channel", "C3", 0.01), mk("channel", "Cz", 0.02), mk("channel", "Pz", 0.30),   # old electrode screen
  mk("channel", "S1_contra", 0.004), mk("channel", "S1_ipsi", 0.03), mk("channel", "SM_furman", 0.2),
  mk("source", "S1", 0.01), mk("source", "ACC", 0.04),
  mk("channel", "S1_contra", NA, "m00_baseline", 0L))
dbWriteTable(con, "gamm_fitted", fitted)

gamm_cfg <- list(channel_units = "clusters", cluster_names = clusters, channel_include = NULL)
h <- suppressMessages(load_hits(con, "v3.1", decomposition = "between_within"))
check(setequal(h$group_value[h$level == "channel"], c("S1_contra", "S1_ipsi")),
      "load_hits (global gamm_cfg = clusters): old electrode hits ignored")
check(setequal(h$group_value[h$level == "source"], c("S1", "ACC")), "load_hits: source hits untouched")
gamm_cfg$channel_units <- "electrodes"
h <- suppressMessages(load_hits(con, "v3.1", decomposition = "between_within"))
check(setequal(h$group_value[h$level == "channel"], c("C3", "Cz")), "load_hits: electrodes")
gamm_cfg$channel_units <- "all"
h <- suppressMessages(load_hits(con, "v3.1", decomposition = "between_within"))
check(setequal(h$group_value[h$level == "channel"], c("C3", "Cz", "S1_contra", "S1_ipsi")), "load_hits: all")
h <- suppressMessages(load_hits(con, "v3.1", decomposition = "between_within", cfg = NULL))
check(nrow(h) == 6, "load_hits(cfg = NULL): every stored hit")
gamm_cfg <- list(channel_units = "clusters", cluster_names = clusters, channel_include = "S1_contra")
h <- suppressMessages(load_hits(con, "v3.1", decomposition = "between_within", levels = "channel"))
check(identical(h$group_value, "S1_contra"), "load_hits: clusters + channel_include")

# -- compare_levels.R path ----
fits <- fitted %>% filter(n_metrics == 1)
f2 <- suppressMessages(filter_channel_units(fits, gamm_cfg, "screened model"))
check(setequal(f2$group_value, c("S1_contra", "S1", "ACC")), "filter_channel_units on screen rows (compare_levels)")
check(nrow(filter_channel_units(fits[0, ], gamm_cfg)) == 0, "filter_channel_units on empty table")
dbDisconnect(con)

# -- analysis_levels: the master switch ----
check(identical(get_analysis_levels(list()), c("channel", "source")), "analysis_levels missing = both (old configs keep working)")
check(identical(get_analysis_levels(list(analysis_levels = "Channel")), "channel"), "analysis_levels case-insensitive")
err <- tryCatch({ get_analysis_levels(list(analysis_levels = "clusters")); FALSE }, error = function(e) TRUE)
check(err, "analysis_levels typo -> error")
ch_only <- list(analysis_levels = "channel")
check(isTRUE(check_screen_level("channel", ch_only)), "screen allowed at a listed level")
err <- tryCatch({ check_screen_level("source", ch_only); FALSE }, error = function(e) grepl("switched off", conditionMessage(e)))
check(err, "source screen refused when analysis_levels = channel")

con <- dbConnect(RSQLite::SQLite(), ":memory:")
dbWriteTable(con, "gamm_fitted", fitted)
gamm_cfg <- list(analysis_levels = "channel", channel_units = "clusters", cluster_names = clusters)
h <- suppressMessages(load_hits(con, "v3.1", decomposition = "between_within"))
check(all(h$level == "channel") && setequal(h$group_value, c("S1_contra", "S1_ipsi")),
      "load_hits default: clusters only, no source hits (holdout, sex diff, figures)")
h <- suppressMessages(load_hits(con, "v3.1", decomposition = "between_within", levels = c("source", "channel")))
check(all(h$level == "channel"), "load_hits explicit levels intersected (robustness$levels)")
h <- suppressMessages(load_hits(con, "v3.1", decomposition = "between_within", levels = "source"))
check(nrow(h) == 0, "load_hits(levels = 'source') with channel-only analysis -> no rows, no error")
gamm_cfg$analysis_levels <- c("channel", "source")
h <- suppressMessages(load_hits(con, "v3.1", decomposition = "between_within"))
check(setequal(h$group_value, c("S1_contra", "S1_ipsi", "S1", "ACC")), "both levels: clusters + source ROIs")
dbDisconnect(con)

# -- classical tests: level and units ----
cfg_c <- list(analysis_levels = "channel", channel_units = "clusters", cluster_names = clusters,
              classical = list(level = NULL, units = c("S1_contra", "S1_ipsi"), rois = c("S1", "ACC")))
check(classical_level(cfg_c) == "channel", "classical level defaults to channel when listed")
check(identical(select_classical_units(avail, "channel", cfg_c), c("S1_contra", "S1_ipsi")), "classical units = S1 clusters")
cfg_c$classical$units <- NULL
check(identical(select_classical_units(avail, "channel", cfg_c), clusters), "classical units NULL = all selected clusters")
cfg_c$classical$level <- "source"
err <- tryCatch({ classical_level(cfg_c); FALSE }, error = function(e) TRUE)
check(err, "classical level not in analysis_levels -> error")
check(identical(select_classical_units(c("ACC", "S1", "M1"), "source", cfg_c), c("S1", "ACC")), "source level uses classical$rois")

cat(sprintf("\n%d / %d checks passed\n", n_pass, n_all))
if (n_pass != n_all) quit(status = 1)