# =============================================================================
# build_95ByBP_epoch_behaviour.R                                      V1.0.0
# -----------------------------------------------------------------------------
# Align 95ByBP behaviour to the EXISTING epochs. No preprocessing rerun.
#
# Why this works: preproc_core.m epochs every laser event (S 32..S 35) after
# ICA and never rejects epochs, so epoch k = laser event k in marker order.
# The 80 rated trials (authors' rule: next kept event is S100) are a subset
# of those epochs. This script writes a behaviour file with ONE ROW PER
# EPOCH, so the merge's positional pairing (trial = epoch index) is correct
# with no change to the merge:
#   - rated, kept by authors    -> laser_power + Painrating from the release
#   - rated, deleted by authors -> Painrating NA (rating not released)
#   - unrated (calibration)     -> Painrating NA
# laser_power is filled for every epoch (marker code S 32..35 -> the
# subject's 1st..4th power level), so intensity checks see no gaps.
#
# Safety: a subject is written only if
#   (1) its epoch count == its laser-event count in the .vmrk files,
#   (2) the epoch event codes match the marker codes one-to-one, and
#   (3) rule C holds (rated markers minus the authors' bad epochs reproduce
#       the released intensity sequence exactly).
# Anything else (e.g. pop_epoch dropped an epoch at a session boundary) is
# left out and reported in the QC file.
#
# Inputs : raw_root (BIDS .vmrk), epoch_root (searched recursively for
#          sub-XXX_epoch_metadata.csv), 95ByBP_behaviour_reindexed.csv
#          (from reindex_95ByBP_behaviour.R).
# Outputs: out_file (goes in behav_dir as 95ByBP_behaviour.csv) and
#          <out_file>_qc.csv (one row per subject).
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

# -- Paths (edit) ----
raw_root    <- "/cifs/seminowicz/eegPainDatasets/CNED/95ByBP"
epoch_root  <- "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/95ByBP"
reidx_file  <- "95ByBP_behaviour_reindexed.csv"
out_file    <- "95ByBP_behaviour_epochs.csv"

laser_codes <- c("S 32" = 1L, "S 33" = 2L, "S 34" = 3L, "S 35" = 4L)
kept_codes  <- c("S  8", "S 12", "S 16", "S 20", "S 32", "S 33", "S 34", "S 35",
                 "S 64", "S 65", "S 66", "S 67", "S100", "S128", "S129", "S130", "S131")
n_rated_expected <- 80L

# -- Markers (same reader as check_95ByBP_events.R V2) ----
read_vmrk <- function(f) {
  L <- grep("^Mk[0-9]+=", readLines(f, warn = FALSE), value = TRUE)
  p <- strsplit(sub("^Mk[0-9]+=", "", L), ",")
  data.frame(type = vapply(p, `[`, "", 1), code = trimws(vapply(p, `[`, "", 2)),
             stringsAsFactors = FALSE)
}
natural_order <- function(x) {
  nums <- regmatches(x, gregexpr("\\d+", x))
  pad  <- vapply(nums, function(n) paste(sprintf("%08d", as.integer(n)), collapse = "_"), "")
  order(pad, x)
}
laser_events <- function(files) {
  ev <- bind_rows(lapply(files[natural_order(files)], read_vmrk)) %>%
    filter(type == "Stimulus", code %in% kept_codes) %>%
    mutate(next_code = lead(code))
  ev %>% filter(code %in% names(laser_codes)) %>%
    mutate(rank = unname(laser_codes[code]),
           rated = !is.na(next_code) & next_code == "S100")
}
norm_code <- function(x) gsub("\\s+", " ", trimws(x))   # "S  32" / "S 32" -> "S 32"

# -- Inputs ----
reidx <- read.csv(reidx_file, stringsAsFactors = FALSE)
vmrk  <- list.files(raw_root, pattern = "\\.vmrk$", recursive = TRUE, full.names = TRUE)
vsubj <- sub(".*(sub-\\d+).*", "\\1", vmrk)
emeta <- list.files(epoch_root, pattern = "^sub-\\d+_epoch_metadata\\.csv$",
                    recursive = TRUE, full.names = TRUE)
esubj <- sub("_epoch_metadata\\.csv$", "", basename(emeta))
cat(length(unique(vsubj)), "subjects with markers,", length(emeta), "epoch metadata files\n")

rows <- list(); qc <- list()
for (s in sort(unique(esubj))) {
  id  <- as.integer(sub("sub-", "", s))
  q   <- data.frame(subj = s, n_epochs = NA, n_laser = NA, n_rated = NA,
                    codes_match = NA, rule_C = NA, status = NA_character_)
  em  <- read.csv(emeta[esubj == s][1], stringsAsFactors = FALSE)
  em  <- em[order(em$epoch_index), ]
  q$n_epochs <- nrow(em)
  if (!any(vsubj == s)) { q$status <- "no_vmrk"; qc[[s]] <- q; next }
  lz <- laser_events(vmrk[vsubj == s])
  q$n_laser <- nrow(lz); q$n_rated <- sum(lz$rated)
  b  <- reidx[reidx$ID == id, ]; b <- b[order(b$Trial_Num), ]
  if (nrow(b) != n_rated_expected) { q$status <- "no_reindexed_behaviour"; qc[[s]] <- q; next }

  # (1)+(2) epochs == laser events, same codes in the same order
  q$codes_match <- nrow(em) == nrow(lz) && all(norm_code(em$event_code) == lz$code)
  if (!q$codes_match) { q$status <- "epochs_differ_from_markers"; qc[[s]] <- q; next }
  if (q$n_rated != n_rated_expected) { q$status <- "rated_count_not_80"; qc[[s]] <- q; next }

  # (3) rule C: rated marker ranks at kept positions == released power ranks
  lv   <- sort(unique(b$laser_power[!is.na(b$laser_power)]))
  keep <- !b$authors_bad_epoch
  rated_rank <- lz$rank[lz$rated]
  q$rule_C <- length(lv) == 4 && identical(rated_rank[keep], match(b$laser_power[keep], lv))
  if (!q$rule_C) { q$status <- "rule_C_failed"; qc[[s]] <- q; next }

  rated_trial <- rep(NA_integer_, nrow(lz)); rated_trial[lz$rated] <- seq_len(n_rated_expected)
  out <- data.frame(ID = id, Trial_Num = em$epoch_index,
                    laser_power = lv[lz$rank],
                    Painrating = NA_real_,
                    event_code = lz$code, rated_trial = rated_trial,
                    authors_bad_epoch = FALSE)
  r <- which(lz$rated)
  out$Painrating[r]        <- b$Painrating
  out$authors_bad_epoch[r] <- b$authors_bad_epoch
  rows[[s]] <- out; q$status <- "ok"; qc[[s]] <- q
}

qc  <- bind_rows(qc)
res <- bind_rows(rows)
write.csv(res, out_file, row.names = FALSE, na = "")
write.csv(qc, sub("\\.csv$", "_qc.csv", out_file), row.names = FALSE)

cat("\nStatus per subject:\n"); print(table(qc$status))
cat("\nRows written:", nrow(res), "epochs for", n_distinct(res$ID), "subjects\n")
cat("  with a rating          :", sum(!is.na(res$Painrating)), "\n")
cat("  rated, authors deleted :", sum(res$authors_bad_epoch), "\n")
cat("  unrated (calibration)  :", sum(is.na(res$rated_trial)), "\n")
cat("Wrote", out_file, "and its _qc.csv\n")