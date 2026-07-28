# =============================================================================
# split_helpers.R
# -----------------------------------------------------------------------------
# Incremental train/holdout split assignment.
#
# Design constraints (agreed):
#   - Stratified by cap_size (32/62/64-channel), proportional split
#     (default 75/25), so every cap size is represented in both train and
#     holdout — required for the channel-pooled GAMM's random effect to be
#     evaluable on held-out subjects (s(channel, bs="re") can't predict for
#     channel levels never seen in training).
#   - Incremental: assigned the moment a subject is first ingested into
#     `subjects`. A subject's split_group is set ONCE and never changed by
#     later runs, even as more subjects of the same cap_size arrive and the
#     running proportions shift.
#   - Reproducible: fixed seed (merge_cfg$split_seed), but because
#     assignment is incremental and order-dependent (each new subject's
#     train/holdout draw depends on how many of their cap_size are already
#     assigned), re-running from scratch with the same seed and the same
#     ingestion order reproduces the same assignments. Changing ingestion
#     order (e.g. processing experiments in a different sequence) WILL
#     change individual assignments on a from-scratch rebuild — this is a
#     real limitation, not a bug, and is why existing assignments are never
#     overwritten: the source of truth is what's already in subject_split,
#     not what a fresh run would compute.
# =============================================================================

library(DBI)
library(dplyr)

# =============================================================================
# init_split_table()
# -----------------------------------------------------------------------------
# Creates subject_split if it doesn't exist. No-op otherwise.
# =============================================================================
init_split_table <- function(con) {
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS subject_split (
      global_subjid VARCHAR PRIMARY KEY,
      cap_size      VARCHAR,
      split_group   VARCHAR,
      split_seed    INTEGER,
      assigned_at   VARCHAR
    )
  ")
}

# =============================================================================
# assign_new_subjects_to_split()
# -----------------------------------------------------------------------------
# Quota-based assignment: for each cap_size group, tracks the running
# train/holdout counts (existing assignments + new ones being made this
# run) and assigns each new subject to whichever group keeps the running
# proportion closest to train_frac. This converges to exactly train_frac
# over time, unlike independent random draws (which can drift several
# percentage points from the target for smaller cap_size groups — verified:
# for n=124 at p=0.75, independent draws had a 5th-95th percentile range of
# roughly 68%-81%, wide enough to matter for the smallest cap_size group's
# holdout stability).
#
# Within each cap_size group, the ORDER subjects are processed in still
# needs to be randomized (seeded) — otherwise the first ceil(n*0.75)
# subjects in global_subjid order always land in train, which isn't a
# random split, just a deterministic one. So: shuffle new subjects within
# each cap_size group (seeded), then walk through in that shuffled order
# assigning each to whichever group is currently furthest below its target
# share.
# =============================================================================
assign_new_subjects_to_split <- function(con, train_frac, seed) {
  init_split_table(con)

  all_subjects <- dbGetQuery(con, "SELECT global_subjid, cap_size FROM subjects")
  existing <- dbGetQuery(con, "SELECT global_subjid, cap_size, split_group FROM subject_split")

  new_subjects <- all_subjects %>% filter(!global_subjid %in% existing$global_subjid)

  if (nrow(new_subjects) == 0) {
    message("  No new subjects to assign — subject_split is up to date.")
    return(invisible(NULL))
  }

  message("  Assigning split for ", nrow(new_subjects), " new subject(s)...")

  set.seed(seed + nrow(existing))

  assigned_rows <- list()

  for (this_cap in unique(new_subjects$cap_size)) {

    # Running counts so far for this cap_size, from EXISTING assignments
    # only (never recomputed from scratch — this is what makes the split
    # stable across incremental runs).
    existing_this_cap <- existing %>% filter(cap_size == this_cap)
    n_train_so_far   <- sum(existing_this_cap$split_group == "train")
    n_holdout_so_far <- sum(existing_this_cap$split_group == "holdout")

    new_this_cap <- new_subjects %>%
      filter(cap_size == this_cap) %>%
      arrange(global_subjid) %>%
      slice_sample(prop = 1)  # seeded shuffle within this cap_size group

    for (i in seq_len(nrow(new_this_cap))) {
      total_so_far <- n_train_so_far + n_holdout_so_far
      current_train_share <- if (total_so_far == 0) 0 else n_train_so_far / total_so_far

      # Assign to train if doing so keeps us at-or-below target share;
      # otherwise assign to holdout. This is a simple greedy quota-fill
      # that converges to train_frac exactly as n grows, without ever
      # revisiting earlier assignments.
      if (current_train_share < train_frac) {
        grp <- "train"
        n_train_so_far <- n_train_so_far + 1
      } else {
        grp <- "holdout"
        n_holdout_so_far <- n_holdout_so_far + 1
      }

      assigned_rows[[length(assigned_rows) + 1]] <- tibble(
        global_subjid = new_this_cap$global_subjid[i],
        cap_size      = this_cap,
        split_group   = grp,
        split_seed    = seed,
        assigned_at   = as.character(Sys.time())
      )
    }
  }

  new_assignments <- bind_rows(assigned_rows)
  dbWriteTable(con, "subject_split", new_assignments, append = TRUE)

  tab <- new_assignments %>% count(cap_size, split_group)
  message("  New assignments by cap_size:")
  for (i in seq_len(nrow(tab))) {
    message("    ", tab$cap_size[i], " / ", tab$split_group[i], ": ", tab$n[i])
  }

  # Report the resulting overall proportion per cap_size (existing + new)
  # so a drift from train_frac is visible immediately, not discovered later.
  all_split <- dbGetQuery(con, "SELECT cap_size, split_group FROM subject_split")
  summary_tab <- all_split %>%
    count(cap_size, split_group) %>%
    group_by(cap_size) %>%
    mutate(share = round(n / sum(n), 3)) %>%
    ungroup()
  message("  Running totals (all assignments to date):")
  for (i in seq_len(nrow(summary_tab))) {
    message("    ", summary_tab$cap_size[i], " / ", summary_tab$split_group[i],
            ": ", summary_tab$n[i], " (", summary_tab$share[i], ")")
  }

  invisible(NULL)
}