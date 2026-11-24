# =============================================================================
# check_95ByBP_events.R                                              V2.0.0
# -----------------------------------------------------------------------------
# Read-only check: which laser events in 95ByBP are the rated behavioural
# trials? Run before changing preprocessing for this experiment.
#
# Background (Nov 2026): every 95ByBP subject failed the merge alignment
# check (99-197 laser epochs vs 80 behavioural rows). The .vmrk markers show
#   - laser stimuli S 32..S 35 = intensity 1..4 (other S-codes are non-laser
#     stimuli of the multi-stimulus paradigm, already not epoched);
#   - S100 = rating marker, ~4-13 s after a rated stimulus;
#   - unrated laser stimuli (runs at one intensity, 15-270 s apart: likely
#     calibration) are followed by another stimulus, not by S100.
# The behaviour file has trailing laser_power 0 / rating 0 rows (padding to
# 80, no stimulus) and NaN ratings (stimulus delivered, no rating).
#
# For each subject this script
#   1. reads ALL its .vmrk files (runs in natural sort order),
#   2. takes the laser events and flags "rated" = next stimulus marker is
#      S100 within max_gap_s,
#   3. maps behaviour laser_power to intensity rank 1..4 (non-zero levels,
#      sorted) after dropping trailing 0/0 padding,
#   4. compares the EEG intensity sequence with the behaviour sequence under
#      two rules:
#        A: rated laser events  == all non-padding behaviour rows
#           (S100 follows even when the rating was not recorded / NaN)
#        B: rated laser events  == behaviour rows with a rating (non-NaN)
#           (no S100 when the participant did not rate)
#        C: rated laser events with the authors' bad epochs removed
#           == all non-padding rows                       (V2, Nov 2026)
#   and reports, per subject, which rule (if any) gives an EXACT match.
#
# V2 (Nov 2026): "rated" now mirrors preprocessing_95.m exactly: keep only
# the authors' event codes (drops S 4, S228-231, ...), then a laser event is
# a trial when the NEXT kept event is S100 (no time limit, sessions merged
# in order). Rule C uses the bad-epoch list from
# all_experiments_preprocessing.csv: the authors deleted those epochs (on the
# 1..80 rated-trial index) after attaching ratings, and the released
# behaviour file is the post-deletion list padded with 0/0 rows.
#
# Exact sequence match over 40-80 trials of 4 intensities is strong evidence
# the pairing is right; a subject matching neither rule stays excluded.
#
# Run (R, anywhere):  source("check_95ByBP_events.R")
# Writes: check_95ByBP_events.csv in the working directory.
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

raw_root   <- "/cifs/seminowicz/eegPainDatasets/CNED/95ByBP"
behav_file <- "/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics/R/data/behaviour/95ByBP_behaviour.csv"
laser_codes <- c("S 32" = 1L, "S 33" = 2L, "S 34" = 3L, "S 35" = 4L)
rating_code <- "S100"
prep_file   <- "all_experiments_preprocessing.csv"   # authors' bad-epoch log
# Event codes the authors keep (pop_selectevent in preprocessing_95.m)
kept_codes  <- c("S  8", "S 12", "S 16", "S 20", "S 32", "S 33", "S 34", "S 35",
                 "S 64", "S 65", "S 66", "S 67", "S100", "S128", "S129", "S130", "S131")

# -- Markers ----
read_vmrk <- function(f) {
  vhdr <- sub("\\.vmrk$", ".vhdr", f)
  si <- grep("^SamplingInterval=", readLines(vhdr, warn = FALSE), value = TRUE)
  fs <- 1e6 / as.numeric(sub(".*=", "", si[1]))
  L  <- grep("^Mk[0-9]+=", readLines(f, warn = FALSE), value = TRUE)
  p  <- strsplit(sub("^Mk[0-9]+=", "", L), ",")
  data.frame(type = vapply(p, `[`, "", 1), code = trimws(vapply(p, `[`, "", 2)),
             t_s = as.numeric(vapply(p, `[`, "", 3)) / fs, file = basename(f),
             stringsAsFactors = FALSE)
}
natural_order <- function(x) {             # ses-2 before ses-10
  nums <- regmatches(x, gregexpr("\\d+", x))
  pad  <- vapply(nums, function(n) paste(sprintf("%08d", as.integer(n)), collapse = "_"), "")
  order(pad, x)
}

rated_laser_sequence <- function(files) {
  # Sessions merged in order, as pop_mergeset does; only the authors' codes
  # survive (pop_selectevent), so S 4 / S228-231 etc. never sit between a
  # laser event and its S100.
  ev <- bind_rows(lapply(files[natural_order(files)], read_vmrk)) %>%
    filter(type == "Stimulus", code %in% kept_codes) %>%
    mutate(next_code = lead(code), next_gap = lead(t_s) - t_s)
  laser <- ev %>% filter(code %in% names(laser_codes)) %>%
    mutate(intensity = unname(laser_codes[code]),
           rated = !is.na(next_code) & next_code == rating_code)
  laser
}

# -- Authors' bad epochs (1..80 rated-trial index) ----
read_bad_epochs <- function(f, block_pattern = "95 By BP") {
  raw <- read.csv(f, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = "character", fileEncoding = "UTF-8-BOM")
  b0  <- grep(block_pattern, unlist(raw[1, ]), fixed = TRUE)
  cols <- b0:(b0 + 3)
  hdr <- tolower(trimws(unlist(raw[2, cols])))
  blk <- raw[-(1:2), cols, drop = FALSE]; names(blk) <- hdr
  blk <- blk[grepl("^\\s*\\d+\\s*$", blk$sub_id), ]
  parse <- function(s) { s <- trimws(s); if (s %in% c("", "\\", ".")) integer(0)
                         else as.integer(regmatches(s, gregexpr("\\d+", s))[[1]]) }
  setNames(lapply(blk[["bad epoch"]], parse), as.integer(blk$sub_id))
}
bad <- if (file.exists(prep_file)) read_bad_epochs(prep_file) else {
  warning("prep_file not found: rule C skipped"); list() }

# -- Behaviour ----
beh <- read.csv(behav_file, stringsAsFactors = FALSE)
behaviour_sequences <- function(b) {
  b <- b[order(b$Trial_Num), ]
  r <- suppressWarnings(as.numeric(b$Painrating))
  pad <- 0
  for (i in rev(seq_len(nrow(b)))) if (b$laser_power[i] == 0 && isTRUE(r[i] == 0)) pad <- pad + 1 else break
  b <- b[seq_len(nrow(b) - pad), ]; r <- r[seq_len(nrow(b))]
  lv <- sort(unique(b$laser_power[b$laser_power > 0]))
  rank <- match(b$laser_power, lv)
  list(all = rank, rated = rank[!is.na(r)], n_pad = pad, n_levels = length(lv))
}

first_mismatch <- function(a, b) {
  n <- min(length(a), length(b)); i <- which(a[seq_len(n)] != b[seq_len(n)])
  if (length(i)) i[1] else if (length(a) != length(b)) n + 1 else NA_integer_
}

# -- Per subject ----
vmrk_all <- list.files(raw_root, pattern = "\\.vmrk$", recursive = TRUE, full.names = TRUE)
subj_of  <- sub(".*(sub-\\d+).*", "\\1", vmrk_all)
cat(length(vmrk_all), ".vmrk files for", n_distinct(subj_of), "subjects\n")

res <- bind_rows(lapply(sort(unique(subj_of)), function(s) {
  id <- as.integer(sub("sub-", "", s))
  b  <- beh[beh$ID == id, ]
  if (nrow(b) == 0) return(data.frame(subj = s, status = "no_behaviour"))
  bs <- behaviour_sequences(b)
  lz <- tryCatch(rated_laser_sequence(vmrk_all[subj_of == s]), error = function(e) NULL)
  if (is.null(lz)) return(data.frame(subj = s, status = "vmrk_read_error"))
  eeg <- lz$intensity[lz$rated]
  be  <- bad[[as.character(id)]]
  eeg_c <- if (is.null(be)) NULL else if (length(be)) eeg[-be[be <= length(eeg)]] else eeg
  data.frame(
    subj = s, n_files = sum(subj_of == s), n_laser = nrow(lz), n_rated_eeg = length(eeg),
    n_beh_rows = length(bs$all), n_beh_rated = length(bs$rated), n_padding = bs$n_pad,
    n_levels = bs$n_levels,
    match_A_all_rows = identical(eeg, bs$all), match_B_rated_rows = identical(eeg, bs$rated),
    n_bad_epochs = if (is.null(be)) NA_integer_ else length(be),
    match_C_bad_removed = !is.null(eeg_c) && identical(eeg_c, bs$all),
    first_mismatch_A = first_mismatch(eeg, bs$all), first_mismatch_B = first_mismatch(eeg, bs$rated),
    first_mismatch_C = if (is.null(eeg_c)) NA_integer_ else first_mismatch(eeg_c, bs$all),
    status = "checked")
}))

write.csv(res, "check_95ByBP_events.csv", row.names = FALSE)
ok <- res$status == "checked"
cat("\nSubjects checked:", sum(ok), "\n")
cat("Exact match, rule A (rated laser events == all non-padding rows):", sum(res$match_A_all_rows[ok]), "\n")
cat("Exact match, rule B (rated laser events == rows with a rating)  :", sum(res$match_B_rated_rows[ok]), "\n")
cat("Exact match, rule C (rated events minus authors' bad epochs)    :", sum(res$match_C_bad_removed[ok]), "\n")
cat("Rated laser events == 80                                        :", sum(res$n_rated_eeg[ok] == 80), "\n")
cat("No rule                                                         :",
    sum(!res$match_A_all_rows[ok] & !res$match_B_rated_rows[ok] & !res$match_C_bad_removed[ok]), "\n")
cat("Other statuses:", paste(capture.output(print(table(res$status[!ok]))), collapse = " "), "\n")
cat("\nFirst 10 subjects:\n")
print(head(res[, c("subj", "n_files", "n_laser", "n_rated_eeg", "n_beh_rows", "n_beh_rated",
                   "n_bad_epochs", "match_A_all_rows", "match_C_bad_removed", "first_mismatch_C")], 10),
      row.names = FALSE)