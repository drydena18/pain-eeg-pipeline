# Preprocessing Stages
This page documents every stage executed by `preproc_core.m`, its inputs, outputs, config parameters, default values, and resume behaviour.

---

## Overview
Stages run in a fixed order. Each stage is independently enabled/disabled in the JSON config. Completed stages are detected by the presence of their output `.set` file and skipped on rerun.

| # | Stage | Dir | Default Tag |
|---|-------|-----|-------------|
– | Import + Montage | *(no save)* | – |
1 | Filter | `01_filter` | `fir` |
2 | Notch | `02_notch` | `notch60` |
3 | Resample | `03_resample` | `rs500` |
4 | Re-Reference | `04_reref` | `reref` |
5 | INITREJ | `05_initrej` | `initrej` |
6 | ICA (+ICLabel) | `06_ica` | `ica` |
7 | Epoch | `07_epoch` | `epoch` |
8 | Baseline Correction | `08_base` | `base` |

---

## Pre-Stage: Import and Montage

Not saved. Runs before the first stage on every execution.

Raw file resolution (`resolve_raw_file.m`)  
Three resolution strategies are tried in order:
1. **BIDS Path** – `P.INPUT.EXP/sub-XXX/eeg/sub-XXX_task-<task>_eeg.bdf`
2. **JSON Pattern** – `cfg.exp.raw.pattern` is `sprintf`-formatted with `subjid` and resolved relative to `P.INPUT.EXP`
3. **Recursive Search** – scans `P.INPUT.EXP/&ast;&ast;` for `&ast;sub-XXX&ast;.bdf/eeg` (if `cfg.exp.raw.search_recursive = true`)

### Import
```matlab
EEG = pop_biosig(rawPath);
EEG = normalize_chan_labels(EEG); % strip spaces/dashes from labels
```

### Montage (optional)
Triggered by `cfg.exp.montage.enabled = true`.
1. Reads a CSV with columns `raw_label`, `std_label`
2. Optionally filters to A1-A32 and B1-B32 channels only (`select_ab_only`, default `true`)
3. Relabels channels; validates no duplicates introduced
4. Optionally applies coordinate lookup from `P.CORE.ELP_FILE` (`do_lookup`, default `true`)
5. Writes audit TSV: `LOGS/sub-XXX_channelmap_applied.tsv`

### Config
```json
"montage": {
    "enabled": true,
    "csv": "montage.csv",
    "select_ab_only": true,
    "do_lookup": true
}
```

### Deterministic RNG
```matlab
rng(double(subjid), 'twister')
```
Set before any stage runs. Ensures ICA produces identical results across reruns and machines for same subject.

---

## Stage 1: Filter
**Function**: `pop_eegfiltnew` (FIR)  
Applies a bandpass filter. Both highpass and lowpass are applied in one call.

### Config
```json
"filter": {
    "enabled": true,
    "tag": "fir",
    "type": "fir",
    "highpass_hz": 0.5,
    "lowpass_hz": 40.0
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
`highpass_hz` | `0.5` | High-pass cutoff (Hz) |
`lowpass_hz` | `40.0` | Low-pass cutoff (Hz) |
`type` | `"fir"` | Filter type label (informational; EEGLAB uses FIR via `pop_eegfiltnew`) |

**Output**: `01_filter/<prefix><subjid>_fir.set`

---

## Stage 2: Notch
**Function**: `pop_eegfiltnew` with `revfilt = 1` (bandstop)  
Removes line noise. Implemented as a bandstop (notch) filter centered at `freq_hz` with half_bandwidth `bw_hz`.

### Config:
```json
"notch": {
    "enabled": true,
    "tag": "notch60",
    "freq_hz": 60,
    "bw_hz": 2
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
`freq_hz` | `60` | Notch center frequency (Hz) |
`bw_hz` | `2` | Half-bandwidth; stop bands is `[freq_hz - bw_hz, freq_hz + bw_hz]`

**Output**: `02_notch/<prefix><subjid>_fir_notch60.set`

---

## Stage 3: Resample
**Function**: `pop_resample`  
Resamples the data to a target sampling rate. If the data is already at `target_hz`, this stage is skipped (no file saved, no tag appended)  
### Config:
```json
"resample": {
    "enabled": false,
    "tag": "rs500",
    "target_hz": 500
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
`target_hz` | `[]` |  Target sampling rate (Hz). Must be set if enabled.

**Output**: `03_resample/<prefix><subjid>_fir_notch60_rs500.set`

---

##  Stage 4: Re-Reference
**Function**: `pop_reref`

Re-references the data. Supports average reference or a specified channel set.
### Config:
```json
"reref": {
    "enabled": true,
    "tag": "reref",
    "mode": "average",
    "channels": []
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
`mode` | `"average"` | `"average"` or `"channels"` |
`channels` | `[]` | Channel indices/labels for `mode = "channels"`. Must be non-empty if mode is `"channels"`

**Output**: `04_reref/<prefix><subjid>_..._reref.set`

---

## Stage 5: INITREJ (Bad channel detection + Interpolation)
The most interactive stage. Combines automated suggestions with a mandatory manual decision.

**Automated Suggestions** (`suggest_bad_channels`)  
Two `pop_rejchan` passes:
- Probability z-score threshold => 5 (normalized)
- Kurtosis z-score threshold => 5 (normalized)

Channels flagged by both passes have their reasons concatenated in the log.

**Channel PSD metrics** (`compute_channel_psd_metrics`)  
Computed for all channels and saved to `LOGS/sub-XXX_chan_psd_metrics.csv`:

| Metric | Description |
|--------|-------------|
`line_ratio` | Power at 59–61 Hz relative to flanking bands (55–59, 61–65 Hz) |
`hf_ratio` | Power 20–40 Hz relative to 1–12 Hz |
`drift_ratio` | Power 1–2 Hz relative to 1–12 Hz |
`alpha_ratio` | Power 8–12 Hz relative to 1–40 Hz |

### QC Plots (saved to LOGS)
- Histogram and bar chart of per-channel STD and RMS
- Topoplot of STD, RMS, and all four PSD metrics
- Median <u>+</u> IQR PSD across all channels
- PSD overlay for suggested bad channels only

**Manual Prompt**:
```text
[INITREJ] Suggested bad channels: [3 17 42]
    1) 3 (Fp1): probability z > 5 (pop_rejchan)
    2) 17 (F7): kurtosis z > 5 (pop_rejchan)
    3) 42 (T7): probability > 5 + kurtosis > 5

Type channel indices to INTERPOLATE (e.g., [1 2 17])
Default = none (Press Enter or type []).
Channels to interpolate:
```
The operator reviews QC plots, then enters the indices to interpolate (or nothing). Suggestions are advisory, the operator has final say. Both the interpolated channels and any suggested-but-not-interpolated channels are logged.

**Config**:
```json
"initrej": {
    "enabled": true,
    "tag": "initrej",
    "badchan": {
        "enabled": false
    },
    "badseg": {
        "enabled": false
    }
}
```

**Output**: `05_initrej/<prefix><subjid>_..._initrej.set`

---

## Stage 6. ICA + ICLabel
**ICA training copy**</br>
Before running ICA, an optional training copy is created with bad segments removed. Bad segments are detected by a single amplitude threshold (`badseg.threshold_uv`) applied across all channels.

The operator is prompted to confirm removal:
```text
Remove detected bad segments from ICA training copy? (y/n) [n]:
```

If confirmed, `pop_select` removes the intervals from the training copy. ICA weights are then transferred back to the full-length EEG.

### ICA
**Function**: `pop_runica`
```json
"ica": {
    "enabled": true,
    "tag": "ica",
    "method": "runica",
    "runica": {
        "extended": false
    }
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
`method` | `"runica"` | ICA algorithm passed to `pop_runica` |
`runica.extended` | `false` | Enabled extended infomax (sub- and super-Gaussian sources) |

### ICLabel (optional)
```json
"iclabel": {
    "enabled": true,
    "tag": "iclabel",
    "auto_reject": false,
    "thresholds": {
        "eye": 0.8,
        "muscle": 0.8,
        "heart": 0.8,
        "line_noise": 0.8,
        "channel_noise": 0.8
    }
}
```

ICLabel classifies each IC into 7 categories: Brain, Muscle, Eye, Heart, Line Noise, Channel Noise, Other. Any IC exceeding a threshold for any non-brain category is flagged for review.

For each flagged IC, a QC packet is saved to `QC/sub-XXX_icqc/`:
- Scalp map (topoplot)
- PSD (0–80 Hz, log scale)
- Activation time series (first 10s)
- ICLabel probabilities in the figure title

Per-IC metrics are also computed and saved to `LOGS/sub-XXX_ic_psd_metrics.csv`.  
The operator is then prompted:
```text
[ICREJ] ICLabel suggested ICs: [1 3 7]
Review QC figs in QC/sub-XXX_icqc/ before deciding.
Type IC indices to REMOVE (e.g., [1 3 7]).
Default = remove none (press Enter or type []).
ICs to remove:
```

**Output**: `06_ica/<prefix><subjid>_..._ica.set` (or `..._ica_iclabel.set` if `iclabel.tag` is set)

---

## Stage 7. Epoch
**Function**: `pop_epoch`

Epochs the continuous data around specified event types.

Before epoching, `validate_events_before_epoch` logs all unique event types present in the data and confirms at least one requested type is found. An error is thrown if none are found.

**Config**:
```json
"epoch": {
    "enabled": true,
    "tag": "epoch",
    "event_types": ["S1", "S2"],
    "tmin_sec": -1.0,
    "tmax_sec": 2.0
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
`event_types` | `{}` | **Required** if enabled. List of event type strings to epoch around. |
`tmin_sec` | `-1.0` | Epoch start relative to event (seconds) |
`tmax_sec` | `2.0` | Epoch end relative to event (seconds) |

**Output**: `07_epoch/<prefix><subjid>_..._epoch.set`

---

## Stage 8: Baseline Correction
**Function**: `pop_rmbase`

Subtracts the mean amplitude within a baseline window from every epoch.

**Config**:
```json
"baseline": {
    "enabled": true,
    "tag": "base",
    "window_sec": [-0.5, 0]
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
`window_sec` | `[0.5, 0]` | Baseline window in seconds. Converted to millisecodns for `pop_rmbase`

**Output**: `08_base/<prefix><subjid>_..._base.set`

---

## Resume Logic
All stages follow the same pattern:
```matlab
[EEG, tags, didLoad] = maybe_load_stage(stageDir, P, subjid, tags, nextTag, logf, EEG);
if ~didLoad
    %...process...
    tags{end+1} = nextTag;
    save_stage(stageDir, P, subjid, tags, EEG, logf);
end
```

`maybe_load_stage` constructs the expected filename (using current tags + `nextTag`) and checks if it exists. If found, it loads the file, updates `tags`, and returns `didLoad = true`. No processing occurs.

This means that if you kill a pipeline mid-run, you can simply rerun it and it will continue from the last completed stage. You can also re-run a single subject by passing `subjects_override = [subjid]`.
> **NOTE**: To force a stage re-run, delete the corresponding `.set` file from the stage directory.

---

## Full Config Schema Reference
```json
{
  "exp": {
    "id": "exp01",
    "out_prefix": "26BB_62_",
    "subjects": [1, 2, 3, 4],
    "raw": {
        "pattern": "sub-%03d.bdf",
        "search_recursive": true
    },

    "montage": {
        "enabled": false,
        "csv": "montage.csv",
        "select_ab_only": true,
        "do_lookup": true
    },

    "channel_locs": {
        "use_elp": false
    },

    "task": ""
  },

  "preproc": {
    "filter": {
        "enabled": true,
        "tag": "fir",
        "type": "fir",
        "highpass_hz": 0.5,
        "lowpass_hz": 40.0
    },

    "notch": {
        "enabled": true,
        "tag": "notch60",
        "freq_hz": 60,
        "bw_hz": 2.0
    },

    "resample": {
        "enabled": false,
        "tag": "rs500",
        "target_hz": 500
    },

    "reref": {
        "enabled": true,
        "tag": "reref",
        "mode": "average",
        "channels": [],
        "badchan": {
            "enabled": false
        },
        "badseg": {
            "enabled": false,
            "threshold_uv": 150
        }
    },

    "ica": {
        "enabled": true,
        "tag": "ica",
        "method": "runica",
        "runica": {
            "extended": false
        },
        "iclabel": {
            "enabled": false,
            "tag": "iclabel",
            "auto_reject": false,
            "thresholds": {
                "eye": 0.8,
                "muscle": 0.8,
                "heart": 0.8,
                "channel": 0.8,
                "line": 0.8
            }
        }
    },

    "epoch": {
        "enabled": true,
        "tag": "epoch",
        "event_types": ["S1"],
        "tmin_sec": -1.0,
        "tmax_sec": 2.0
    },

    "baseline": {
        "enabled": true,
        "tag": "base",
        "window_sec": [-0.5, 0]
    }
  }
```