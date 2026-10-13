# =============================================================================
# alignment_helpers.R
# -----------------------------------------------------------------------------
# Behaviour <-> EEG trial alignment check (merge stage).
#
# Behaviour is joined to EEG features on (global_subjid, trial), where the
# EEG "trial" is the epoch index. That is only correct if the k-th epoch is
# the k-th behavioural trial. Any extra or missing epoch (a mis-specified
# event code, a trial lost at a session boundary, a missing session) shifts
# every later rating onto the wrong trial without any error.
#
# assess_alignment() returns one row per subject:
#   status = "ok"                 counts agree (and intensity check passed or
#                                 was not applicable)
#            "count_mismatch"     epochs != behavioural trials, or the
#                                 spectral / source files disagree with the
#                                 epoch count
#            "intensity_mismatch" low/high epoch labels disagree with laser_power
#            "behaviour_noncontiguous"  behavioural trials are not 1..n
#            "no_eeg"             behaviour only (nothing processed yet)
#            "no_behaviour"       EEG only
#   intensity_check = "pass" | "fail" | "labels_unknown" | "single_level" |
#                     "no_metadata" | "not_run"
# =============================================================================

library(dplyr)

assess_alignment <- function(behav, epochs, n_spec, n_src) {
  out <- list(n_behav = NA_integer_, n_epochs = NA_integer_,
              intensity_check = "not_run", intensity_agreement = NA_real_,
              status = "ok", reason = "")

  has_behav <- !is.null(behav) && nrow(behav) > 0
  has_meta  <- !is.null(epochs) && nrow(epochs) > 0
  has_eeg   <- has_meta || !is.na(n_spec) || !is.na(n_src)

  if (has_behav) out$n_behav <- nrow(behav)
  if (has_meta)  out$n_epochs <- nrow(epochs)

  if (!has_eeg) {
    out$status <- "no_eeg"; out$reason <- "no EEG outputs found"
    return(out)
  }
  if (!has_behav) {
    out$status <- "no_behaviour"; out$reason <- "EEG outputs but no behavioural rows"
    return(out)
  }
  if (!identical(sort(behav$trial), seq_len(nrow(behav)))) {
    out$status <- "behaviour_noncontiguous"
    out$reason <- "behavioural trial numbers are not 1..n"
    return(out)
  }

  # Reference EEG trial count: epoch metadata if present, else feature files
  n_ref <- if (has_meta) out$n_epochs else if (!is.na(n_spec)) n_spec else n_src
  ref_label <- if (has_meta) "epochs" else if (!is.na(n_spec)) "spectral" else "source"

  bad <- character(0)
  if (!is.na(n_spec) && n_spec != n_ref) bad <- c(bad, sprintf("spectral %d", n_spec))
  if (!is.na(n_src)  && n_src  != n_ref) bad <- c(bad, sprintf("source %d", n_src))
  if (length(bad) > 0) {
    out$status <- "count_mismatch"
    out$reason <- sprintf("%s has %d trials but %s", ref_label, n_ref, paste(bad, collapse = ", "))
    return(out)
  }
  if (n_ref != out$n_behav) {
    out$status <- "count_mismatch"
    out$reason <- sprintf("%d EEG trials (%s) vs %d behavioural trials", n_ref, ref_label, out$n_behav)
    return(out)
  }

  # Intensity check (only meaningful once counts agree)
  if (!has_meta) {
    out$intensity_check <- "no_metadata"
    return(out)
  }
  lab <- tolower(epochs$intensity)
  if (!all(c("low", "high") %in% lab)) {
    out$intensity_check <- if (all(lab %in% c("unknown", "", NA))) "labels_unknown" else "single_level"
    return(out)
  }
  j <- inner_join(epochs %>% select(trial, intensity),
                  behav  %>% select(trial, laser_power), by = "trial") %>%
    filter(tolower(intensity) %in% c("low", "high"), !is.na(laser_power))
  lv <- sort(unique(j$laser_power))
  if (length(lv) != 2) {
    out$intensity_check <- "not_run"
    out$reason <- sprintf("low/high labels present but %d laser levels", length(lv))
    return(out)
  }
  expected <- ifelse(tolower(j$intensity) == "high", lv[2], lv[1])
  out$intensity_agreement <- mean(j$laser_power == expected)
  if (out$intensity_agreement == 1) {
    out$intensity_check <- "pass"
  } else {
    out$intensity_check <- "fail"
    out$status <- "intensity_mismatch"
    out$reason <- sprintf("low/high labels match laser_power on %.0f%% of trials",
                          100 * out$intensity_agreement)
  }
  out
}