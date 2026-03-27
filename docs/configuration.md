# Configuration Reference

This page is the single authoritative reference for the JSON configuration files that drive the pipeline. Every field is listed with its type, default value, and a description of its effect.

Config files live in `config/` and are named by experiment: `config/exp01.json`.

---

## How Configuration Works

The JSON file is loaded by `load_cfg()` and passed as a struct through the call chain:

```
expXX.json → load_cfg() → preproc_default() → normalize_preproc_defaults() → preproc_core()
```

`preproc_default` validates required fields and resolves subjects. `normalize_preproc_defaults` fills every missing field with its default. By the time `preproc_core` runs, the config struct is fully populated regardless of how minimal the JSON was.

Fields you omit from the JSON take their default values. Fields you include override the defaults.

---

## Top-Level Structure

```json
{
  "exp":     { ... },
  "preproc": { ... },
  "spectral": { ... },
  "paths":   { ... }
}
```

`"exp"` and `"preproc"` are required. `"spectral"` is required for the spectral pipeline. `"paths"` is optional.

---

## `paths` — Path Overrides

All paths have hardcoded defaults pointing to the lab server. Override here to run on a different machine without editing any MATLAB code.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `raw_root` | string | `/cifs/seminowicz/eegPainDatasets/CNED` | Root directory containing all raw experiment folders |
| `proj_root` | string | `/cifs/seminowicz/eegPainDatasets/CNED/da-analysis` | Root directory for all outputs. Also accepts `out_root` as an alias. |

```json
"paths": {
    "raw_root":  "/Volumes/seminowicz/CNED",
    "proj_root": "/Users/me/analysis"
}
```

---

## `exp` — Experiment Settings

### Required Fields

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Experiment identifier. Must match a key in `experiment_registry()` (e.g. `"exp01"`). Underscores are stripped before lookup. |
| `out_prefix` | string | Filename prefix for all output `.set` files (e.g. `"26BB_62_"`). Must be non-empty. |

### Subject Resolution

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `subjects` | array of int | *(from TSV)* | List of subject IDs to process. If omitted or empty, resolved from `participants.tsv`. |

Subject resolution priority: `subjects_override` argument → `cfg.exp.subjects` → `participants.tsv` in `P.INPUT.EXP` → `participants.tsv` in `P.RESOURCE`.

Subject IDs can be integers or strings like `"sub-001"` — trailing digits are extracted automatically.

### Raw File Resolution

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `raw.pattern` | string | `""` | `sprintf` pattern for the raw filename, formatted with `subjid`. E.g. `"sub-%03d_task-pain_eeg.bdf"`. Resolved relative to `P.INPUT.EXP`. |
| `raw.search_recursive` | bool | `true` | If pattern fails, recursively search `P.INPUT.EXP/**` for any `.bdf`/`.eeg` file matching the subject number. |
| `task` | string | `""` | Task label used in BIDS path construction (`sub-XXX_task-<task>_eeg.bdf`). Falls back to `raw_dirname` from the registry if empty. |

### Montage

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `montage.enabled` | bool | `false` | Whether to apply a channel relabelling map |
| `montage.csv` | string | `""` | Path to the montage CSV file (columns: `raw_label`, `std_label`). Relative paths are resolved from `P.RESOURCE`. |
| `montage.select_ab_only` | bool | `true` | If true, keep only channels matching `A1–A32` and `B1–B32` before relabelling (Biosemi 64-channel layout) |
| `montage.do_lookup` | bool | `true` | After relabelling, look up 3D coordinates from `P.CORE.ELP_FILE` via `pop_chanedit` |

### Channel Locations

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `channel_locs.use_elp` | bool | `false` | Load channel coordinates from `P.CORE.ELP_FILE` directly (without montage relabelling). Use when labels already match the ELP file. |

---

## `preproc` — Preprocessing Stages

Every stage block has at minimum an `enabled` bool and a `tag` string. The tag is appended to the output filename when the stage completes.

### `preproc.filter`

FIR bandpass filter via `pop_eegfiltnew`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | |
| `tag` | string | `"fir"` | |
| `type` | string | `"fir"` | Informational label only; EEGLAB always uses FIR via `pop_eegfiltnew` |
| `highpass_hz` | number | `0.5` | High-pass cutoff (Hz) |
| `lowpass_hz` | number | `40` | Low-pass cutoff (Hz) |

---

### `preproc.notch`

Bandstop (notch) filter via `pop_eegfiltnew` with `revfilt = 1`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | |
| `tag` | string | `"notch60"` | |
| `freq_hz` | number | `60` | Notch centre frequency (Hz) |
| `bw_hz` | number | `2` | Half-bandwidth; stop band is `[freq_hz − bw_hz, freq_hz + bw_hz]` |

---

### `preproc.resample`

Resample via `pop_resample`. Skipped (no tag appended) if data is already at `target_hz`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `false` | |
| `tag` | string | `"rs500"` | |
| `target_hz` | number | `[]` | Target sampling rate. **Must be set if enabled.** |

---

### `preproc.reref`

Re-reference via `pop_reref`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | |
| `tag` | string | `"reref"` | |
| `mode` | string | `"average"` | `"average"` for average reference, `"channels"` for specific channels |
| `channels` | array | `[]` | Channel indices or labels for `mode = "channels"`. Error if empty when mode is `"channels"`. |

---

### `preproc.initrej`

Bad channel detection, QC plotting, and manual spherical interpolation.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | |
| `tag` | string | `"initrej"` | |
| `badchan.enabled` | bool | `false` | Reserved; not currently used to gate any automated removal |
| `badseg.enabled` | bool | `false` | Enable bad segment detection for the ICA training copy |
| `badseg.threshold_uv` | number | — | Amplitude threshold (µV) for bad segment detection. Required if `badseg.enabled = true`. |

---

### `preproc.ica`

Independent component analysis via `pop_runica`, optionally followed by ICLabel.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | |
| `tag` | string | `"ica"` | |
| `method` | string | `"runica"` | ICA algorithm; passed directly to `pop_runica` as `'icatype'` |
| `runica.extended` | bool | `false` | Enable extended infomax (models super-Gaussian sources, better for some artifact types) |

#### `preproc.ica.iclabel`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `false` | Run ICLabel classification and generate QC packets |
| `tag` | string | `""` | If non-empty, appended to the output filename after the ICA tag |
| `auto_reject` | bool | `false` | Reserved; not currently implemented. Rejection is always manual. |
| `thresholds.eye` | number | — | Probability threshold for Eye category |
| `thresholds.muscle` | number | — | Probability threshold for Muscle category |
| `thresholds.heart` | number | — | Probability threshold for Heart category |
| `thresholds.line_noise` | number | — | Probability threshold for Line Noise category |
| `thresholds.channel_noise` | number | — | Probability threshold for Channel Noise category |

Any `thresholds` field that is omitted simply disables that category's check. You can enable ICLabel but only check for `eye` and `muscle` by setting only those two.

---

### `preproc.epoch`

Epoch continuous data around event markers via `pop_epoch`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | |
| `tag` | string | `"epoch"` | |
| `event_types` | array of string | `[]` | **Required when enabled.** Event type strings to epoch around (e.g. `["S1", "S2"]`). Error if empty. |
| `tmin_sec` | number | `-1.0` | Epoch start relative to event (seconds) |
| `tmax_sec` | number | `2.0` | Epoch end relative to event (seconds) |

---

### `preproc.baseline`

Baseline correction via `pop_rmbase`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | |
| `tag` | string | `"base"` | |
| `window_sec` | array [lo, hi] | `[-0.5, 0]` | Baseline window in seconds. Converted to milliseconds internally. |

---

## `spectral` — Spectral Pipeline Settings

### Top-Level

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `true` | If false, `spectral_default` exits without running `spectral_core` |
| `input_stage` | string | `"08_base"` | Stage directory to load the input `.set` file from |

### `spectral.psd`

Welch PSD parameters.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `fmin_hz` | number | `1` | Lower frequency bound of PSD output |
| `fmax_hz` | number | `80` | Upper frequency bound of PSD output |
| `window_sec` | number | `2.0` | Welch window length (seconds). Adapted down if shorter than trial length. |
| `overlap_frac` | number | `0.5` | Fractional overlap between Welch windows |
| `nfft` | int | `0` | FFT size. `0` = auto (next power of 2 ≥ window length). |

### `spectral.alpha`

Alpha frequency band definitions.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `alpha_hz` | array [lo, hi] | `[8, 12]` | Total alpha band |
| `slow_hz` | array [lo, hi] | `[8, 10]` | Slow alpha sub-band |
| `fast_hz` | array [lo, hi] | `[10, 12]` | Fast alpha sub-band |

### `spectral.qc`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `plot_mode` | string | `"summary"` | `"summary"` (default), `"debug"` (subset of trials), or `"exhaustive"` (all trials) |
| `save_heatmaps` | bool | `true` | Save the channel × trial heatmap panel figure |
| `legend_max_channels` | int | `20` | Maximum channels shown in plot legends |
| `max_debug_trials` | int | `5` | Number of trials selected for debug mode plots |

### `spectral.fooof`

Python FOOOF bridge. Requires `fooof` / `specparam` installed in the Python environment.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `false` | |
| `python_exe` | string | `"python3"` | Python executable name or full path |
| `script_path` | string | `""` | **Required when enabled.** Absolute path to `fooof_bridge.py` |
| `fmin_hz` | number | `2` | Lower frequency bound for FOOOF fitting |
| `fmax_hz` | number | `40` | Upper frequency bound for FOOOF fitting |
| `aperiodic_mode` | string | `"fixed"` | `"fixed"` or `"knee"` |
| `peak_width_limits` | array [lo, hi] | `[1, 12]` | Min/max peak width (Hz) |
| `max_n_peaks` | int | `6` | Maximum number of peaks fitted per trial |
| `min_peak_height` | number | `0.0` | Minimum peak height above aperiodic component |
| `peak_threshold` | number | `1.0` | Number of SDs above noise floor for peak detection |
| `alpha_band_hz` | array [lo, hi] | `[8, 12]` | Band to search for the alpha peak in FOOOF output |
| `verbose` | bool | `false` | Pass verbose flag to the Python FOOOF object |

### `spectral.trial_spectral`

Pre vs. post-stimulus PSD figures per trial.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `enabled` | bool | `false` | |
| `pre_sec` | array [lo, hi] | `[-1.0, 0.0]` | Pre-stimulus time window (seconds relative to event) |
| `post_sec` | array [lo, hi] | `[0.0, 1.0]` | Post-stimulus time window |
| `max_trials` | int | `20` | Maximum number of trials to save figures for |
| `fmin_hz` | number | `1` | Lower frequency bound of the per-trial PSD |
| `fmax_hz` | number | `80` | Upper frequency bound |
| `legend_max_channels` | int | `16` | Maximum channels shown in the figure legend |

---

## Complete Annotated Example

```json
{
    "paths": {
        "raw_root":  "/cifs/seminowicz/eegPainDatasets/CNED",
        "proj_root": "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis"
    },

    "exp": {
        "id":         "exp01",
        "out_prefix": "26BB_62_",
        "subjects":   [1, 2, 3, 4, 5],

        "raw": {
            "pattern":          "",
            "search_recursive": true
        },

        "montage": {
            "enabled":        true,
            "csv":            "montage.csv",
            "select_ab_only": true,
            "do_lookup":      true
        },

        "channel_locs": {
            "use_elp": false
        }
    },

    "preproc": {
        "filter": {
            "enabled":      true,
            "tag":          "fir",
            "type":         "fir",
            "highpass_hz":  0.5,
            "lowpass_hz":   40
        },

        "notch": {
            "enabled":  true,
            "tag":      "notch60",
            "freq_hz":  60,
            "bw_hz":    2
        },

        "resample": {
            "enabled":    false,
            "tag":        "rs500",
            "target_hz":  500
        },

        "reref": {
            "enabled":  true,
            "tag":      "reref",
            "mode":     "average",
            "channels": []
        },

        "initrej": {
            "enabled": true,
            "tag":     "initrej",
            "badchan": { "enabled": false },
            "badseg":  {
                "enabled":       true,
                "threshold_uv":  150
            }
        },

        "ica": {
            "enabled": true,
            "tag":     "ica",
            "method":  "runica",
            "runica":  { "extended": false },
            "iclabel": {
                "enabled":     true,
                "tag":         "iclabel",
                "auto_reject": false,
                "thresholds": {
                    "eye":           0.80,
                    "muscle":        0.80,
                    "heart":         0.80,
                    "line_noise":    0.80,
                    "channel_noise": 0.80
                }
            }
        },

        "epoch": {
            "enabled":      true,
            "tag":          "epoch",
            "event_types":  ["S1"],
            "tmin_sec":     -1.0,
            "tmax_sec":     2.0
        },

        "baseline": {
            "enabled":     true,
            "tag":         "base",
            "window_sec":  [-0.5, 0]
        }
    },

    "spectral": {
        "enabled":      true,
        "input_stage":  "08_base",

        "psd": {
            "fmin_hz":      1,
            "fmax_hz":      80,
            "window_sec":   2.0,
            "overlap_frac": 0.5,
            "nfft":         0
        },

        "alpha": {
            "alpha_hz": [8, 12],
            "slow_hz":  [8, 10],
            "fast_hz":  [10, 12]
        },

        "qc": {
            "plot_mode":           "summary",
            "save_heatmaps":       true,
            "legend_max_channels": 20,
            "max_debug_trials":    5
        },

        "fooof": {
            "enabled":            true,
            "python_exe":         "python3",
            "script_path":        "/home/user/pain-eeg-pipeline/python/fooof_bridge.py",
            "fmin_hz":            2,
            "fmax_hz":            40,
            "aperiodic_mode":     "fixed",
            "peak_width_limits":  [1, 12],
            "max_n_peaks":        6,
            "min_peak_height":    0.0,
            "peak_threshold":     1.0,
            "alpha_band_hz":      [8, 12],
            "verbose":            false
        },

        "trial_spectral": {
            "enabled":             false,
            "pre_sec":             [-1.0, 0.0],
            "post_sec":            [0.0,  1.0],
            "max_trials":          20,
            "fmin_hz":             1,
            "fmax_hz":             80,
            "legend_max_channels": 16
        }
    }
}
```

---

## Minimal Config (Preprocessing Only)

The smallest valid config that will run a complete preprocessing pipeline. All unlisted fields take their defaults.

```json
{
    "exp": {
        "id":         "exp01",
        "out_prefix": "26BB_62_"
    },
    "preproc": {
        "epoch": {
            "event_types": ["S1"]
        }
    }
}
```

`event_types` must always be specified because there is no safe default for it.