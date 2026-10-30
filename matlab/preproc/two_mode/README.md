# Two-mode preprocessing (`preproc_2mode`)

Runs the same preprocessing as `expXX_preproc` -> `preproc_core`, split so the
slow part runs unattended and the manual IC review happens later, in one sitting.

```matlab
preproc_2mode('exp06')               % unattended: Mode 1 -> Mode 2a (+ Mode 2b for anyone already reviewed)
preproc_2mode('exp06', 'review')     % IC rejection at your desk, then Mode 2b
preproc_2mode('exp06', 'status')     % manifest only
preproc_2mode('exp06', 'run', [1 4]) % subset
```

Re-run any time: finished work is detected from commit markers and skipped.
If MATLAB crashes on participant 17, participants 1-16 are not redone.

## Modes

| Mode | Runs | What it does | Writes |
|---|---|---|---|
| 1 | unattended | 00 concat/load -> 01 filter -> 02 notch -> 03 resample -> 04 reref -> 05 initrej -> clean copy -> rank-aware extended Infomax -> weights transferred | `00_concat`-`05_initrej` (same files as preproc_core), `06_ica_ready/` |
| 2a | unattended | ICLabel, suggested rejections, QC packets, IC PSD metrics, topo grid | `QC/sub-XXX_icqc/`, `LOGS/`, `06_ica_review/` |
| review | at desk | component scroll -> `prompt_ic_reject` -> confirmation scroll | `06_ica_review/*_icdecisions.json` |
| 2b | unattended | remove ICs -> 07 epoch -> 08 baseline -> 09 Hilbert | `06_ica/`, `07_epoch/`, `08_base/`, `09_hilbert/` (same names as preproc_core) |

Spectral (`input_stage: "08_base"`) and source localization read Mode 2b's
output without changes.

## What differs from preproc_core

* **Bad channels (05)** - `prompt_channel_interp` is replaced by
  `interp_policy` in `pp_hooks.m` (`'suggested'` = interpolate every channel
  `suggest_bad_channels` flags; `'none'`). The channel scroll is skipped.
* **Clean copy** - the `prompt_yesno` in `make_ica_training_copy` is replaced
  by `badseg_policy` (`'remove'` | `'keep'`). Same >threshold detection.
* **ICA rank** - passed explicitly (`'pca', r`). Because reref (04) runs before
  interpolation (05): rank = nChan - nInterp, or nChan - 1 if nothing was
  interpolated. Cross-checked against an eigenvalue estimate; the smaller wins.
* **IC review** - after ICA instead of inside it; decisions saved per
  participant and bound to the decomposition (weights hash).
* **Seed** - `rng(<subject number>)` right before `runica`, as preproc_core.

Both policies are part of the Mode 1 config hash and are recorded in
`EEG.etc.initrej`, `EEG.etc.pipeline` and the ICA-ready marker.

## Resume and staleness

* Mode 1 reuses existing `00`-`05` stage files (like preproc_core). A reused
  `05_initrej` from a manual preproc_core run doesn't record which channels
  were interpolated, so the eigenvalue rank is used for that participant.
* Editing the experiment JSON marks the affected mode out of date:
  filter/notch/resample/reref/initrej/ICA -> Mode 1 *stale* (reported, not
  re-run, so reviewed ICs are kept; use `'rerun_stale'`); ICLabel thresholds /
  grid -> Mode 2a re-runs; epoch/baseline/Hilbert -> Mode 2b re-runs.
* Code changes inside `st_*.m` are invisible to the hash: bump
  `pp.modeX.version` in `pp_hooks.m`.

Manifest and errors: `<RUN_ROOT>/_pp_logs/pp_manifest.csv`, `pp_errors.log`.
Per-participant logs: `sub-XXX/LOGS/sub-XXX_pp_<mode>.log`.

## Files

```
matlab/preproc/preproc_2mode.m      entry point (next to expXX_preproc.m)
matlab/preproc/two_mode/            runner: pp_run, pp_review, pp_status, pp_mode1/2a/2b, ...
matlab/preproc/two_mode/pp_hooks.m  wiring + the two unattended policies
matlab/preproc/two_mode/stages/     st_*.m - one per preproc_core block
```
