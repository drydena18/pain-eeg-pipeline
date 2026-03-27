# Helper Function Reference

This page documents every standalone MATLAB helper function in the preprocessing framework. These live as individual `.m` files and are callable from anywhere on the MATLAB path.

Functions are grouped by responsibility. For spectral-specific helpers (`spec_*.m`), see the [Spectral README](spectral-README.md).

> **Note:** `preproc_helpers.m` is not used. All helpers are standalone files.

---

## Filesystem

### `ensure_dir(d)`
Creates directory `d` if it does not exist. Silently no-ops if it already exists. Accepts both `char` and `string` input.

```matlab
ensure_dir(LOGS);
ensure_dir(fullfile(subRoot, '01_filter'));
```

---

### `open_log(LOGS, subjid, stem)` → `[fid, logPath, cleanupObj]`
Opens a per-subject log file for writing. Returns the file descriptor, the resolved path, and an `onCleanup` object that closes the file when the caller goes out of scope.

| Argument | Default | Description |
|----------|---------|-------------|
| `LOGS` | — | Directory path |
| `subjid` | — | Integer subject ID |
| `stem` | `'log'` | Filename stem; produces `sub-XXX_<stem>.log` |

If the file cannot be opened, falls back to stdout (`fid = 1`) with a warning.

---

### `safeClose(fid)`
Closes a file descriptor only if it is not stdout (`fid ~= 1`) and is positive. Used in `onCleanup` handlers throughout the pipeline to avoid closing stdout accidentally.

---

### `save_stage(stageDir, P, subjid, tags, EEG, logf)`
Saves the current EEG struct to a stage directory as a `.set`/`.fdt` pair. The filename is constructed from the current `tags` cell array via `P.NAMING.fname`. Logs the output path.

---

### `maybe_load_stage(stageDir, P, subjid, tags, nextTag, logf, EEG)` → `[EEG, tags, didLoad]`
Checks whether a cached `.set` file exists for the next stage. If found, loads it and returns `didLoad = true`. If not found, returns the input `EEG` and `didLoad = false` so the caller proceeds with processing.

The expected filename is constructed from `tags` + `nextTag`. This is the core of the stage resumption mechanism.

| Argument | Description |
|----------|-------------|
| `stageDir` | Stage output directory |
| `P` | Paths struct from `config_paths` |
| `subjid` | Integer subject ID |
| `tags` | Current cumulative tag cell array |
| `nextTag` | Tag string for the stage being checked |
| `logf` | Log file descriptor |
| `EEG` | Current EEG struct (returned unchanged if cache hit) |

---

### `local_fname(subjid, tags, prefix, defaultPrefix)` → `fname`
Constructs a `.set` filename from its components.

```
<prefix><subjid_zero_padded_3>_<tag1>_<tag2>_..._<tagN>.set
```

`prefix` overrides `defaultPrefix` if non-empty. Both `string` and `char` tags are accepted.

```matlab
fname = local_fname(1, {'fir', 'notch60', 'reref'}, [], '26BB_62_');
% → '26BB_62_001_fir_notch60_reref.set'
```

---

## Logging

### `logmsg(fid, fmt, varargin)`
Writes a timestamped message to a log file and to stdout simultaneously.

```
[YYYY-MM-DD HH:MM:SS] message text
```

`fmt` and `varargin` are passed to `sprintf`. If `fid == 1`, only writes once (avoids double stdout output).

---

## Config and Path Resolution

### `load_cfg(json_path)` → `cfg`
Reads a JSON file and returns a MATLAB struct via `jsondecode`. Throws a clear error if the file does not exist.

---

### `get_cfg_string(cfg, keyPath, defaultVal)` → `val`
Safely traverses a nested struct using a field-name array and returns the value as a `string`. Returns `defaultVal` if any level of the path is missing or empty. Never errors.

```matlab
root = get_cfg_string(cfg, ["paths", "raw_root"], "/default/path");
```

---

### `resolve_raw_file(P, cfg, subjid)` → `rawPath`
Resolves the raw EEG file path for a given subject using three strategies in order:

1. **BIDS path:** `P.INPUT.EXP/sub-XXX/eeg/sub-XXX_task-<task>_eeg.bdf` (and `.BDF`, `.eeg`, `.EEG` variants)
2. **JSON pattern:** `cfg.exp.raw.pattern` is `sprintf`-formatted with `subjid` and resolved relative to `P.INPUT.EXP`
3. **Recursive search:** `dir` glob over `P.INPUT.EXP/**` for `*sub-XXX*.bdf/.eeg`

Returns an empty string if no file is found. The caller in `preproc_core` checks and skips the subject with a warning.

---

### `resolve_exp_out(P, cfg)` → `expOut`
Returns the experiment output directory name. Preference order: `P.EXP.out_dirname` → `cfg.exp.out_dirname` → `cfg.exp.id` → `"experiment"`.

---

### `experiment_registry()` → `R`
Returns a struct of all registered experiments. Each field is an experiment entry with `id`, `raw_dirname`, `out_dirname`, and `out_prefix`. Called by `config_paths`.

---

## Config Normalization

### `normalize_preproc_defaults(Pp)` → `Pp`
Fills all missing fields in `cfg.preproc` with sensible defaults. Called by `preproc_default`. Uses `ensureBlock` and `defaultField` internally.

Default values applied:

| Block | Field | Default |
|-------|-------|---------|
| `filter` | `enabled` | `true` |
| `filter` | `highpass_hz` | `0.5` |
| `filter` | `lowpass_hz` | `40` |
| `notch` | `enabled` | `true` |
| `notch` | `freq_hz` | `60` |
| `notch` | `bw_hz` | `2` |
| `resample` | `enabled` | `false` |
| `reref` | `enabled` | `true` |
| `reref` | `mode` | `"average"` |
| `initrej` | `enabled` | `true` |
| `ica` | `enabled` | `true` |
| `ica` | `method` | `"runica"` |
| `epoch` | `enabled` | `true` |
| `epoch` | `tmin_sec` | `-1.0` |
| `epoch` | `tmax_sec` | `2.0` |
| `baseline` | `enabled` | `true` |
| `baseline` | `window_sec` | `[-0.5, 0]` |

Also validates that `epoch.event_types` is non-empty when `epoch.enabled = true`.

---

### `defaultField(s, field, val)` → `s`
Sets `s.(field) = val` only if the field is missing or empty. Non-destructive — never overwrites an existing value.

```matlab
cfg.filter = defaultField(cfg.filter, 'highpass_hz', 0.5);
```

---

### `defaultStruct(s, field)` → `s`
Sets `s.(field) = struct()` only if the field is missing or empty. Used to ensure sub-structs exist before `defaultField` calls on their children.

---

### `ensureBlock(S, field, defaultEnabled, defaultTag)` → `S`
Ensures `S.(field)` exists as a struct, then applies `defaultField` for both `enabled` and `tag`.

```matlab
Pp = ensureBlock(Pp, 'filter', true, 'fir');
% → Pp.filter.enabled = true  (if missing)
% → Pp.filter.tag = 'fir'     (if missing)
```

---

### `mustHave(s, field, msg)`
Throws an error with `msg` if `s` does not have `field`. Used for required-field validation in `preproc_default` and `spectral_default`.

```matlab
mustHave(cfg, 'exp', 'Missing cfg.exp in JSON.');
```

---

## Subject ID Handling

### `normalize_subject_ids(subs_raw)` → `subs`
Converts a raw subject list — which may be a numeric array, a string array, or a cell array of strings like `"sub-001"` — into a sorted column vector of positive integers. Extracts trailing digits from string IDs, removes duplicates and non-positive values.

```matlab
subs = normalize_subject_ids(["sub-001", "sub-003", "sub-002"]);
% → [1; 2; 3]
```

Throws a descriptive error if any entry cannot be parsed.

---

### `extract_subject_column(T)` → `subs_raw`
Extracts the subject ID column from a `participants.tsv` table. Accepts column names: `subjid`, `participant_id`, or `subject` (case-insensitive). Throws a clear error if none are found.

---

## EEG Struct Utilities

### `normalize_chan_labels(EEG)` → `EEG`
Strips whitespace and dashes from all channel labels in `EEG.chanlocs`. Applied immediately after raw import to normalize Biosemi-style labels like `"A 1"` or `"A-1"` to `"A1"`.

---

### `ensure_etc_path(EEG)` → `EEG`
Ensures `EEG.etc` exists as a struct. Called before writing any metadata into `EEG.etc` to avoid struct-field-on-empty errors.

---

### `has_chanlocs(EEG)` → `tf`
Returns `true` if `EEG.chanlocs` is non-empty and the first entry has an `X` coordinate field. Used to guard `topoplot` calls — topoplots require 3D coordinates from a prior ELP lookup.

---

### `safe_chan_label(EEG, ch)` → `lbl`
Returns the label string for channel index `ch`. Returns `"ChN"` (e.g. `"Ch42"`) if `chanlocs` is missing, the index is out of range, or the label is empty. Never errors.

---

### `vec2str(v)` → `s`
Converts a numeric vector to a compact bracket-notation string, suitable for log messages.

```matlab
vec2str([1 3 17])   % → '[1 3 17]'
vec2str([])         % → '[]'
```

---

## Montage

### `apply_montage_biosemi_from_csv(P, cfg, EEG, logf, subjid, LOGS)` → `EEG`
Applies a channel relabelling map from a CSV file (columns: `raw_label`, `std_label`) to an EEG struct. Full workflow:

1. Resolves the CSV path relative to `P.RESOURCE` if the path is relative
2. Optionally selects only A1–A32 and B1–B32 channels (`select_ab_only`); uses case-insensitive regex matching
3. Reads the CSV; validates that required columns and uniqueness of both label columns
4. Validates all `raw_label` entries are present in the EEG (case-insensitive); missing channels → error
5. Relabels channels in-place; validates no duplicates introduced (case-insensitive)
6. Optionally looks up 3D coordinates from `P.CORE.ELP_FILE` via `pop_chanedit` (`do_lookup`)
7. Writes `LOGS/sub-XXX_channelmap_applied.tsv` as a permanent audit record

---

### `write_channelmap_tsv(LOGS, subjid, EEG)`
Writes a TSV audit file mapping channel index to label. Called by `apply_montage_biosemi_from_csv`. Format: `index \t label`.

---

## Bad Channel Detection

### `suggest_bad_channels(EEG)` → `[badChans, reasons, metrics]`
Automated bad channel suggestion using two `pop_rejchan` passes (probability z-score > 5, kurtosis z-score > 5). Both passes are wrapped in `try/catch` so a plugin failure does not abort the pipeline.

Returns:
- `badChans` — sorted integer vector of flagged channel indices
- `reasons` — cell array of reason strings (one per channel; concatenated if flagged by both passes)
- `metrics` — struct with `chan_rms` `[nChan × 1]` and `chan_std` `[nChan × 1]`

---

## INITREJ Plots

### `make_initrej_plots(LOGS, subjid, EEG, metrics, badChans)`
Orchestrates all INITREJ QC figure generation. Calls `save_hist`, `save_bar`, `save_topo_metric`, `save_channel_psd_overview`, and `save_channel_psd_badchans`. Also writes `sub-XXX_chalabels.txt`.

---

### `save_hist(LOGS, subjid, v, name)`
Saves a histogram of metric vector `v` as `sub-XXX_initrej_hist_<name>.png`.

---

### `save_bar(LOGS, subjid, v, badChans, name)`
Saves a bar chart of metric vector `v` with vertical dashed lines at `badChans` channel indices as `sub-XXX_initrej_bar_<name>.png`.

---

### `save_topo_metric(LOGS, subjid, EEG, v, label)`
Saves a topoplot of per-channel metric `v`. Requires `has_chanlocs(EEG)` to have 3D coordinates. Output: `sub-XXX_initrej_topo_<label_lowercase>.png`.

---

### `save_channel_psd_overview(LOGS, subjid, EEG)`
Computes a Welch PSD for every channel (2-second windows, 50% overlap) and saves a single figure showing the median ± IQR (25th/75th percentile) across channels (0–80 Hz). Output: `sub-XXX_initrej_psd_overview.png`.

---

### `save_channel_psd_badchans(LOGS, subjid, EEG, badChans)`
Saves a PSD overlay plot for the suggested bad channels only. Each channel's PSD is a separate line. Output: `sub-XXX_initrej_psd_badchans.png`.

---

## Channel PSD Metrics

### `compute_channel_psd_metrics(EEG)` → `chanPSD`
Computes four spectral quality metrics per channel via Welch PSD (2-second windows, 50% overlap, FFT size = next power of 2 ≥ window length):

| Field | Band numerator | Band denominator | Interpretation |
|-------|---------------|-----------------|----------------|
| `line_ratio` | 59–61 Hz | 55–59 + 61–65 Hz | Line noise contamination |
| `hf_ratio` | 20–40 Hz | 1–12 Hz | Muscle artifact |
| `drift_ratio` | 1–2 Hz | 1–12 Hz | Slow drift |
| `alpha_ratio` | 8–12 Hz | 1–40 Hz | Alpha band prominence |

Returns a struct with four `[nChan × 1]` vectors.

---

### `write_channel_psd_csv(LOGS, subjid, EEG, chanPSD)`
Writes `sub-XXX_chan_psd_metrics.csv` with columns: `chan_idx`, `label`, `line_ratio`, `hf_ratio`, `drift_ratio`, `alpha_ratio`.

---

### `save_chan_psd_topos(LOGS, subjid, EEG, chanPSD)`
Saves topoplots for all four channel PSD metrics. Guards on `has_chanlocs`. Calls `save_topo_metric` four times.

---

### `bp_psd(f, pxx, band)` → `p`
Trapezoidal integration (`trapz`) of PSD vector `pxx` over frequency vector `f` within `band = [flo, fhi]`. Returns `0` if no frequency bins fall within the band. Used internally by `compute_channel_psd_metrics`.

### `bandpower_from_psd(f, pxx, band)` → `p`
Identical implementation to `bp_psd`. Used by `compute_ic_psd_metrics` and `spec_compute_alpha_features_from_psd`. The two names exist for historical reasons; both are kept to avoid breaking call sites.

---

## ICA Utilities

### `make_ica_training_copy(EEG, cfg, logf)` → `[EEGtrain, segInfo]`
Creates a copy of EEG for ICA training. If `cfg.preproc.initrej.badseg.enabled = true`, detects samples exceeding `badseg.threshold_uv` across any channel, converts them to intervals via `mask_to_intervals`, and prompts the operator to confirm removal via `prompt_yesno`. If confirmed, removes the intervals from the copy via `pop_select`.

ICA weights trained on this copy are transferred back to the full-length EEG in `preproc_core` after training completes.

`segInfo` records: `removed` (bool), `n_intervals`, `pct_time`, `intervals` (`[N × 2]`).

---

### `mask_to_intervals(mask)` → `intervals`
Converts a logical sample mask to an `[N × 2]` matrix of `[start, end]` sample indices for each contiguous `true` run.

```matlab
mask = [0 0 1 1 1 0 1 1 0];
intervals = mask_to_intervals(mask);
% → [3 5; 7 8]
```

---

## IC QC

### `save_ic_qc_packets(QC, subjid, EEG, icList)`
Generates a three-panel QC figure for each IC in `icList`:
- **Top left:** Scalp map (`topoplot` of `EEG.icawinv(:, ic)`)
- **Top right:** Welch PSD of IC activation (0–80 Hz, dB scale)
- **Bottom:** IC activation time series (first 10 seconds)

ICLabel probabilities are shown in the figure title if `EEG.etc.ic_classification.ICLabel.classifications` is available. Output: `QC/sub-XXX_icqc/sub-XXX_icNNN_qc.png`.

---

### `iclabel_suggest_reject(EEG, thr)` → `[suggestICs, reasons]`
Reads `EEG.etc.ic_classification.ICLabel.classifications` and flags any IC whose probability for a non-brain category meets or exceeds the corresponding threshold in `thr`.

Recognized `thr` fields: `eye`, `muscle`, `heart`, `line_noise`, `channel_noise`.

Returns empty arrays if ICLabel output is not present in `EEG.etc`. Safe to call even if the ICLabel plugin was not run.

---

### `compute_ic_psd_metrics(EEG, icList)` → `icMetrics`
Computes per-IC spectral metrics for each IC in `icList`. IC activations are computed as `(EEG.icaweights * EEG.icasphere) * EEG.data`.

Returns a struct array (one entry per IC) with fields:

| Field | Description |
|-------|-------------|
| `ic` | IC index |
| `peak_hz` | Frequency of peak power (0.5–40 Hz) |
| `bp` | Struct: `delta` (1–4), `theta` (4–8), `alpha` (8–12), `beta` (13–30), `gamma` (30–45) Hz band powers |
| `hf_ratio` | Power 20–40 Hz / power 1–12 Hz |
| `line_ratio` | Power 59–61 Hz / (55–59 + 61–65 Hz) |

---

### `write_ic_metrics_csv(LOGS, subjid, icMetrics)`
Writes `sub-XXX_ic_psd_metrics.csv` with columns: `ic`, `peak_hz`, `delta`, `theta`, `alpha`, `beta`, `gamma`, `hf_ratio`, `line_ratio`.

---

### `log_ic_metrics(logf, icMetrics)`
Writes a one-line-per-IC summary of all PSD metrics to the log file.

---

## Manual Prompts

### `prompt_channel_interp(EEG, suggested)` → `interpChans`
Prints suggested bad channels with their labels and waits for operator input. Returns a sorted vector of channel indices to interpolate, or `[]` if Enter is pressed. Input is parsed via `str2num` and accepts standard MATLAB vector syntax (`[1 2 17]`).

---

### `prompt_ic_reject(suggestICs)` → `removedICs`
Prints the ICLabel-suggested ICs and waits for operator input. Returns a sorted vector of IC indices to remove via `pop_subcomp`, or `[]` (the default) if Enter is pressed.

---

### `prompt_yesno(prompt, defaultTF)` → `tf`
Presents a yes/no prompt and returns a logical. Accepts `y`, `yes`, `n`, `no` (case-insensitive). Returns `defaultTF` for empty input or unrecognized responses.

---

## Event Validation

### `validate_events_before_epoch(EEG, wanted, logf)`
Logs all unique event types present in `EEG.event` with their trial counts. Confirms at least one of the `wanted` event type strings is present. Throws an error if none are found, listing both what was requested and what was found.

---

## Spectral Helpers (`spec_*.m`)

Spectral helpers follow the same conventions but are prefixed `spec_` to avoid name collisions with preprocessing helpers. They are documented in the [Spectral README](spectral-README.md). A brief listing:

| Function | Purpose |
|----------|---------|
| `spec_compute_psd_trials` | Welch PSD per channel per trial → `[nChan × nFreq × nTrials]` |
| `spec_compute_alpha_features_from_psd` | Alpha features from PSD tensor |
| `spec_compute_psd_epoch_window` | PSD for one trial and one time window |
| `spec_alpha_bandwidth_proxy` | Weighted SD bandwidth of alpha power per trial |
| `spec_squeeze_ga_features` | Normalize GA feature struct fields to `[nTrials × 1]` |
| `spec_find_latest_set` | Find newest `.set` file for a subject in a stage directory |
| `spec_get_chanlabels` | Extract channel label strings from `EEG.chanlocs` |
| `spec_run_fooof_python` | MATLAB → Python FOOOF bridge |
| `spec_fill_fooof_alpha` | Fill missing FOOOF alpha peaks with PSD-based fallbacks |
| `spec_write_chan_trial_csv` | Write per-channel × per-trial spectral CSV |
| `spec_write_ga_trial_csv` | Write GA × per-trial spectral + FOOOF CSV |
| `spec_plot_summary` | One-page GA PSD + PAF + alpha interaction summary figure |
| `spec_plot_fooof_summary` | Six-panel FOOOF trial-wise summary figure |
| `spec_plot_heatmap_panel` | Channel × trial heatmap panel for all features |
| `spec_plot_heatmap` | Single channel × trial heatmap |
| `spec_plot_debug_trials` | Per-trial multi-channel PSD overlays |
| `spec_plot_trial_spectral_qc` | Pre vs post-stim PSD per trial |
| `spec_plot_trial_prepost_psd` | Render one pre/post PSD figure |
| `spec_select_trials_to_plot` | Auto-select trials for debug plots (extreme SF and PAF) |
| `spec_ensure_dir` | Same as `ensure_dir`, spectral-namespaced |
| `spec_open_log` | Same as `open_log`, spectral-namespaced |
| `spec_safe_close` | Same as `safeClose`, spectral-namespaced |
| `spec_logmsg` | Same as `logmsg`, spectral-namespaced |