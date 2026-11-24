# Test drop_lateral_units(): run from R/tests/ (RSQLite stands in for DuckDB)
library(DBI)
source("../merge/helpers/lateral_helpers.R")
ok <- 0; fail <- 0
chk <- function(cond, msg) { if (isTRUE(cond)) { ok <<- ok + 1 } else { fail <<- fail + 1; cat("FAIL:", msg, "\n") } }

con <- dbConnect(RSQLite::SQLite(), ":memory:")
units <- c("C3", "C4", "CP1", "S1_contra", "S1_ipsi", "SM_furman")
mk <- function(g) expand.grid(global_subjid = g, trial = 1:3, channel = units,
                              metric_name = c("a", "b"), stringsAsFactors = FALSE)
tm <- rbind(mk("E06_S001"), mk("E01_S001")); tm$metric_value <- 1
dbWriteTable(con, "spectral_metrics", tm)
sm <- unique(tm[, c("global_subjid", "channel", "metric_name")]); sm$metric_value <- 2
dbWriteTable(con, "spectral_subject_metrics", sm)

dbBegin(con)
n <- drop_lateral_units(con, "E06_S001")
dbCommit(con)
chk(n == 2 * 3 * 2 + 2 * 2, "deleted rows = trial rows (2 units x 3 trials x 2 metrics) + subject rows (2 x 2)")
left <- dbGetQuery(con, "SELECT DISTINCT channel FROM spectral_metrics WHERE global_subjid = 'E06_S001'")$channel
chk(setequal(left, c("C3", "C4", "CP1", "SM_furman")), "bilateral subject keeps electrodes + SM_furman")
other <- dbGetQuery(con, "SELECT DISTINCT channel FROM spectral_metrics WHERE global_subjid = 'E01_S001'")$channel
chk(setequal(other, units), "other subjects untouched")
chk(drop_lateral_units(con, "E06_S001") == 0, "second call is a no-op")
chk(drop_lateral_units(con, "E99_S999") == 0, "unknown subject: nothing deleted")

# rollback: a failure after the drop leaves the subject's rows intact
dbWriteTable(con, "spectral_metrics", tm, overwrite = TRUE)
dbBegin(con); invisible(drop_lateral_units(con, "E06_S001")); dbRollback(con)
chk(nrow(dbGetQuery(con, "SELECT * FROM spectral_metrics WHERE channel = 'S1_contra' AND global_subjid = 'E06_S001'")) == 6,
    "rolled-back transaction restores lateral rows")
dbDisconnect(con)
cat(sprintf("\n%d passed, %d failed\n", ok, fail)); if (fail > 0) quit(status = 1)