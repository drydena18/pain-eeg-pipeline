# =============================================================================
# reindex_95ByBP_behaviour.R                                          V1.0.0
# -----------------------------------------------------------------------------
# Rebuild the 95ByBP behaviour file on the ORIGINAL 80-trial index.
#
# What the released file is (Nov 2026 finding):
#   The authors (preprocessing_95.m) merge ses-1..6, keep only the listed
#   stimulus codes, label a laser event L1..L4 when the next kept event is
#   S100, epoch -> 80 trials, attach ratings BY POSITION, then delete bad
#   trials by hand. The released 95ByBP_behaviour.csv is the post-deletion
#   trial list, padded at the end with laser_power 0 / rating 0 rows.
#   The bad-epoch indices in all_experiments_preprocessing.csv are on the
#   pre-deletion 1..80 index: for all 95 subjects, number of trailing 0/0
#   rows == number of listed bad epochs.
#
# What this script does, per subject:
#   original_trial = sort(setdiff(1:80, bad_epochs)) assigned to the kept rows
#   in order; the deleted trials are inserted back as rows with NA
#   laser_power / Painrating and authors_bad_epoch = TRUE.
#   Output has exactly 80 rows per subject, Trial_Num = original trial index,
#   so it pairs with the k-th rated laser epoch (rule C in
#   check_95ByBP_events.R).
#
# Ratings travel with their epochs (attached before deletion), so kept rows'
# ratings are unchanged; only their trial numbers move.
#
# Run:  Rscript reindex_95ByBP_behaviour.R <behaviour.csv> <preprocessing.csv> <out.csv>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
behav_file <- if (length(args) >= 1) args[1] else "/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics/R/data/behaviour/95ByBP_behaviour.csv"
prep_file  <- if (length(args) >= 2) args[2] else "/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics/R/data/behaviour/all experiments_preprocessing.csv"
out_file   <- if (length(args) >= 3) args[3] else "95ByBP_behaviour_reindexed.csv"
n_trials   <- 80L
`%||%` <- function(a, b) if (is.null(a)) b else a

# -- Bad epochs for 95ByBP from the wide preprocessing log --------------------
# Layout: 4 columns per experiment block (sub_id, bad epoch, bad electrode,
# bad ica; column order can differ by block, so match on the header).
read_bad_epochs <- function(f, block_pattern = "95 By BP") {
  raw <- read.csv(f, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = "character", fileEncoding = "UTF-8-BOM")
  b0  <- grep(block_pattern, unlist(raw[1, ]), fixed = TRUE)
  if (length(b0) != 1) stop("Block '", block_pattern, "' not found in row 1 of ", f)
  cols <- b0:(b0 + 3)
  hdr  <- tolower(trimws(unlist(raw[2, cols])))
  blk  <- raw[-(1:2), cols, drop = FALSE]; names(blk) <- hdr
  blk  <- blk[grepl("^\\s*\\d+\\s*$", blk$sub_id), ]
  parse <- function(s) { s <- trimws(s); if (s %in% c("", "\\", ".")) integer(0)
                         else as.integer(regmatches(s, gregexpr("\\d+", s))[[1]]) }
  setNames(lapply(blk[["bad epoch"]], parse), as.integer(blk$sub_id))
}

bad  <- read_bad_epochs(prep_file)
beh  <- read.csv(behav_file, stringsAsFactors = FALSE)
beh$Painrating <- suppressWarnings(as.numeric(beh$Painrating))

out <- list(); log <- list()
for (id in sort(unique(beh$ID))) {
  b <- beh[beh$ID == id, ]; b <- b[order(b$Trial_Num), ]
  # trailing padding: laser_power 0 and rating 0
  pad <- 0L
  for (i in rev(seq_len(nrow(b))))
    if (b$laser_power[i] == 0 && isTRUE(b$Painrating[i] == 0)) pad <- pad + 1L else break
  kept <- b[seq_len(nrow(b) - pad), ]
  be   <- bad[[as.character(id)]] %||% integer(0)
  ok   <- length(be) == pad && nrow(b) == n_trials &&
          all(be >= 1 & be <= n_trials) && !anyDuplicated(be)
  log[[length(log) + 1]] <- data.frame(ID = id, n_rows = nrow(b), n_padding = pad,
                                       n_bad_epochs = length(be), consistent = ok)
  if (!ok) next
  orig <- setdiff(seq_len(n_trials), be)
  full <- data.frame(ID = id, Trial_Num = seq_len(n_trials),
                     laser_power = NA_real_, Painrating = NA_real_,
                     authors_bad_epoch = seq_len(n_trials) %in% be,
                     released_row = NA_integer_)
  full$laser_power[orig]  <- kept$laser_power
  full$Painrating[orig]   <- kept$Painrating
  full$released_row[orig] <- kept$Trial_Num
  out[[length(out) + 1]] <- full
}

log <- do.call(rbind, log)
res <- do.call(rbind, out)
write.csv(res, out_file, row.names = FALSE, na = "")
cat("Subjects:", nrow(log), "| consistent (padding == bad epochs):", sum(log$consistent), "\n")
if (any(!log$consistent)) { cat("Inconsistent, NOT written:\n"); print(log[!log$consistent, ], row.names = FALSE) }
cat("Rows written:", nrow(res), "(", nrow(res) / n_trials, "subjects x", n_trials, ")\n")
cat("Deleted-by-authors trials re-inserted:", sum(res$authors_bad_epoch), "\n")
cat("Kept trials with NaN rating (stimulus delivered, no rating stored):",
    sum(!res$authors_bad_epoch & is.na(res$Painrating)), "\n")
cat("Wrote", out_file, "\n")