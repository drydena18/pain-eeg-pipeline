# Spectral Pipeline

This page documents the trial-wise spectral feature extraction pipeline: its call chain, the features it computes, the FOOOF bridge, and its outputs.

---

## Overview

The spectral pipeline takes baseline-corrected epoched EEG (the output of stage 08) and extracts trial-wise spectral features for every channel and for the grand-average (GA) across channels. It optionally runs FOOOF/specparam via a Python bridge to decompose the GA PSD into aperiodic and periodic (peak) components.

---

## Call Chain

```
exp01_spectral.m
    └── spectral_default.m      validate + normalize config, resolve subjects
          └── spectral_core.m   trial-wise PSD, features, FOOOF, CSV, plots
```

### `exp01_spectral.m`

Experiment entrypoint. Adds EEGLAB and project helpers to the path, loads the JSON config, builds `P` via `config_paths()`, and dispatches to `spectral_default()`.

### `spectral_default.m`

Validates required fields, resolves subjects (same logic as `preproc_default`), fills all defaults for `cfg.spectral.*`, and dispatches to `spectral_core()`. If `cfg.spectral.enabled = false`, it exits without running the core.

### `spectral_core.m`

The main execution loop. For each subject:

1. Resolve and load the input `.set` file from the specified stage directory
2. Compute per-channel × per-trial PSD (`spec_compute_psd_trials`)
3. Compute alpha features from PSD (`spec_compute_alpha_features_from_psd`)
4. Average across channels for grand-average (GA) features
5. Optionally run FOOOF on GA PSD via Python bridge
6. Optionally fill missing FOOOF alpha peak values with PSD-based fallbacks
7. Write per-channel × per-trial CSV and GA × per-trial CSV
8. Generate QC figures

---

## Input

The pipeline reads from `cfg.spectral.input_stage` (default: `"08_base"`).

It searches the stage directory for the newest `.set` file matching the subject prefix and subject ID. The file must be epoched (`EEG.trials > 1`); continuous data will be skipped with a warning.

---

## PSD Computation

**Function:** `spec_compute_psd_trials`

Welch PSD computed per channel, per trial using `pwelch`.

| Config Parameter | Default | Description |
|-----------------|---------|-------------|
| `psd.window_sec` | `2.0` | Welch window length (seconds) |
| `psd.overlap_frac` | `0.5` | Fractional overlap between windows |
| `psd.nfft` | `0` (auto) | FFT size; if 0, uses next power of 2 ≥ window length |
| `psd.fmin_hz` | `1` | Lower frequency bound of output |
| `psd.fmax_hz` | `80` | Upper frequency bound of output |

**Output shape:** `[nChan × nFreq × nTrials]`

The window is adapted to the trial length if the trial is shorter than `window_sec`.

---

## Alpha Features

**Function:** `spec_compute_alpha_features_from_psd`

Computed from the PSD tensor `[nChan × nFreq × nTrials]` across three overlapping frequency bands:

| Config | Default | Band |
|--------|---------|------|
| `alpha.alpha_hz` | `[8, 12]` | Total alpha |
| `alpha.slow_hz` | `[8, 10]` | Slow alpha |
| `alpha.fast_hz` | `[10, 12]` | Fast alpha |

**Features computed (all `[nChan × nTrials]`):**

| Feature | Description |
|---------|-------------|
| `paf_cog_hz` | Peak Alpha Frequency — centre of gravity (CoG) of power within the alpha band |
| `pow_slow_alpha` | Absolute slow alpha power (trapz integration) |
| `pow_fast_alpha` | Absolute fast alpha power |
| `pow_alpha_total` | Total alpha power |
| `rel_slow_alpha` | Slow alpha / total alpha |
| `rel_fast_alpha` | Fast alpha / total alpha |
| `sf_ratio` | Slow / fast power ratio |
| `sf_logratio` | log(slow) − log(fast) |
| `sf_balance` | (slow − fast) / (slow + fast) — bounded [−1, +1] |
| `slow_alpha_frac` | slow / (slow + fast) — bounded [0, 1] |

All power values use trapezoidal integration (`trapz`) over the frequency axis. PAF CoG is computed as the power-weighted mean frequency within the alpha band.

---

## Grand Average (GA) Features

The GA PSD is computed by averaging the `[nChan × nFreq × nTrials]` tensor across channels (power domain, `mean`). This produces `[nFreq × nTrials]`.

The same `spec_compute_alpha_features_from_psd` function is applied to the GA PSD, producing scalar features `[nTrials × 1]` after squeezing via `spec_squeeze_ga_features`.

---

## FOOOF Bridge

**Config:**
```json
"fooof": {
    "enabled": false,
    "python_exe": "python3",
    "script_path": "/path/to/fooof_bridge.py",
    "fmin_hz": 2,
    "fmax_hz": 40,
    "aperiodic_mode": "fixed",
    "peak_width_limits": [1, 12],
    "max_n_peaks": 6,
    "min_peak_height": 0.0,
    "peak_threshold": 1.0,
    "alpha_band_hz": [8, 12],
    "verbose": false
}
```

When enabled, `spec_run_fooof_python` writes three files to `SPECTRAL/tmp/`:

- `sub-XXX_fooof_freqs.csv` — frequency vector
- `sub-XXX_fooof_psd.csv` — GA PSD per trial (`[nTrials × nFreq]`)
- `sub-XXX_fooof_cfg.json` — FOOOF parameters

It then calls the Python script via `system()`:
```
python3 fooof_bridge.py --freq ... --psd ... --cfg ... --out ...
```

The Python script fits a FOOOF model to each trial independently and writes `sub-XXX_fooof_out.json`.

### Python script: `fooof_bridge.py`

Uses the `fooof` (specparam) library. For each trial:

1. Fits the FOOOF model to the GA PSD within `[fmin_hz, fmax_hz]`
2. Extracts aperiodic parameters: offset, exponent
3. Finds the strongest peak within `alpha_band_hz`: CF, power, bandwidth
4. Writes per-trial results as JSON; failed trials produce NaN entries with a `fail_reason` string

### FOOOF Output Filling (`spec_fill_fooof_alpha`)

When FOOOF fails to detect an alpha peak (returns NaN for CF), a fallback is applied:

1. **PSD peak-pick** — finds the frequency of maximum power within the alpha band; estimates bandwidth via FWHM
2. **PAF CoG fallback** — if the PSD window is too short, uses the CoG from `spec_compute_alpha_features_from_psd`

Each filled value carries a `_source` string: `"fooof"`, `"psd_peak"`, or `"paf_cog"`. This provenance information is written to the output CSV.

---

## Output CSVs

### `csv/sub-XXX_spectral_chan_by_trial.csv`

One row per `(trial, channel)` combination.

| Column | Description |
|--------|-------------|
| `subjid` | Subject ID |
| `trial` | Trial number |
| `chan_idx` | Channel index |
| `chan_label` | Channel label |
| `paf_cog_hz` | PAF CoG for this channel × trial |
| `pow_slow_alpha` | Slow alpha power |
| `pow_fast_alpha` | Fast alpha power |
| `pow_alpha_total` | Total alpha power |
| `rel_slow_alpha` | Relative slow alpha |
| `rel_fast_alpha` | Relative fast alpha |
| `sf_ratio` | Slow/fast ratio |
| `sf_logratio` | Slow/fast log ratio |
| `sf_balance` | Slow/fast balance |
| `slow_alpha_frac` | Slow alpha fraction |

### `csv/sub-XXX_spectral_ga_by_trial.csv`

One row per trial (GA across channels).

Includes all GA spectral features plus FOOOF outputs:

| Column | Description |
|--------|-------------|
| `fooof_offset` | Aperiodic offset |
| `fooof_exponent` | Aperiodic exponent (1/f slope) |
| `fooof_r2` | Model R² |
| `fooof_error` | Model fit error |
| `fooof_alpha_cf` | Alpha peak CF from FOOOF (NaN if no peak found) |
| `fooof_alpha_pw` | Alpha peak power from FOOOF |
| `fooof_alpha_bw` | Alpha peak bandwidth from FOOOF |
| `fooof_alpha_cf_filled` | CF after fallback filling |
| `fooof_alpha_cf_source` | Provenance: `"fooof"`, `"psd_peak"`, or `"paf_cog"` |
| `fooof_alpha_pw_filled` | Power after fallback filling |
| `fooof_alpha_pw_source` | Provenance |
| `fooof_alpha_bw_filled` | Bandwidth after fallback filling |
| `fooof_alpha_bw_source` | Provenance |
| `fooof_alpha_found` | 1 if any alpha peak was identified |

---

## Output Figures

### `figures/sub-XXX_spectral_summary.png`

Four-panel summary figure:
- Trial-averaged GA PSD with alpha/slow/fast band markers
- PAF (CoG) over trials
- `sf_balance` and `sf_logratio` over trials
- *(FOOOF summary saved separately as `sub-XXX_spectral_fooof.png`)*

### `figures/sub-XXX_spectral_fooof.png`

Six-panel FOOOF trial-wise summary:
- Aperiodic exponent over trials
- Aperiodic offset over trials
- Alpha peak CF over trials
- Alpha peak power over trials
- Alpha peak bandwidth over trials
- R² and error over trials

### `figures/sub-XXX_heatmap_panel.png`

Channel × trial heatmaps for all spectral features. Channels on the Y-axis, trials on the X-axis. One panel per feature. Useful for identifying trial-by-trial and channel-by-channel patterns.

### `trial_spectral/sub-XXX_trial-NNN_prepost_psd.png`

(Requires `cfg.spectral.trial_spectral.enabled = true`)

Pre-stimulus vs. post-stimulus PSD for each trial, channels overlaid. Median PSD shown in bold. Log Y-axis.

---

## Plot Modes

| Mode | Description |
|------|-------------|
| `"summary"` | Default. Only summary and heatmap figures. |
| `"debug"` | Summary + per-trial PSD channel overlays for a subset of selected trials (extreme sf_balance and PAF values). |
| `"exhaustive"` | Summary + per-trial PSD channel overlays for every trial. |

---

## Full Spectral Config Schema

```json
"spectral": {
    "enabled": true,
    "input_stage": "08_base",

    "psd": {
        "fmin_hz": 1,
        "fmax_hz": 80,
        "window_sec": 2.0,
        "overlap_frac": 0.5,
        "nfft": 0
    },

    "alpha": {
        "alpha_hz": [8, 12],
        "slow_hz":  [8, 10],
        "fast_hz":  [10, 12]
    },

    "qc": {
        "plot_mode": "summary",
        "save_heatmaps": true,
        "legend_max_channels": 20,
        "max_debug_trials": 5
    },

    "fooof": {
        "enabled": false,
        "python_exe": "python3",
        "script_path": "",
        "fmin_hz": 2,
        "fmax_hz": 40,
        "aperiodic_mode": "fixed",
        "peak_width_limits": [1, 12],
        "max_n_peaks": 6,
        "min_peak_height": 0.0,
        "peak_threshold": 1.0,
        "alpha_band_hz": [8, 12],
        "verbose": false
    },

    "trial_spectral": {
        "enabled": false,
        "pre_sec":  [-1.0, 0.0],
        "post_sec": [0.0,  1.0],
        "max_trials": 20,
        "fmin_hz": 1,
        "fmax_hz": 80,
        "legend_max_channels": 16
    }
}
```