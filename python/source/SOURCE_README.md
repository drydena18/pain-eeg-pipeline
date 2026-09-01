# CNED EEG Source Localization Pipeline (Python/MNE)
**Author:** Dryden Arseneau  
**Affiliations:** Schabrun Lab · Seminowicz Lab  
**Dataset:** CNED – Zhao et al., Sci Data (2025)  
https://doi.org/10.1030/s41597-025-05900-1  

---

## 1. Overview
This pipeline performs **sLORETA (standardized Low Resolution Electromagnetic Tomography) source localization** on preprocessed EEG epochs and extracts a comprehensive set of whole-epoch, pre-stimulus, post-stimulus, and LEP (laser-evoked potential) metrics in source space. Field names match the MATLAB channel-space pipeline exactly (see `spectral_core.m`), so it feeds directly into the R GAMM and classical test pipelines for a like-for-like channel-vs-source comparison, joined by metric name.

Because CNED participants have no subject MRI, all localization uses the **MNI152 fsaverage template brain** with a pre-computed BEM (boundary element model) head model. This is standard practice for scalp EEG pain research and introduces ~1-2 cm spatial uncertainty that averages out at the group level.

The pipeline is the Python counterpart to the MATLAB preprocessing and spectral pipelines. It follows the same design principles: config-driven, modular, logged, resumable, and human-auditable.

---

## 2. Design Principles
- **Separation of concerns** – each `src_*.py` module does exactly one thing
- **Compute/write separation** – `src_compute_*` functions return data; `src_write_*` functions write files; neither does both
- **Graceful per-subject failure** – one bad subject never kills the batch; errors are logged and skipped
- **Forward solution caching** – expensive per-subject forward models are saved to disk and reloaded on re-runs
- **Config-driven defaults** – all parameters have safe defaults; only `fsaverage.subjects_dir`, `exp.out_prefix`, and `exp.subjects` must be set explicitly
- **Flat module layout** – all `src_*.py` files live alongside `source_core.py`; no `helpers/` package structure needed

---

## 3. Call Chain

```
expXX_source.py     ← Experiment entry point (CLI + importable)
    ↓
source_default.py   ← Config validation, defaults, run header, dispatch
    ↓
source_core.py      ← Per-subject orchestration loop
    ↓
src_io.py           ← Log open/write/close, .set fle resolution, epoch loading
src_assets.py       ← fsaverage BEM/trans/src validation, label + ROI building
src_inverse.py      ← Forward solution, noise covariance, inverse operator, apply
src_spectral.py     ← Welch PSD, bandpower, CoG, bandpass filter primitives
src_alpha_features.py ← Shared whole/pre/post alpha feature computation (11 metrics/window)
src_erd.py          ← ERD family (fractional pre->post change) + p5_flag QC gate
src_compute_metric_deltas.py ← Generic delta_<name> = post - pre, all 11 metric families
src_add_prefix.py   ← Prefix dict/row-list keys with whole_/pre_/post_
src_merge_rows.py   ← Merge multiple row-lists on (trial, roi_idx) or (roi_idx)
src_prestim.py      ← Pre-stimulus Hilbert phase (t=0) + TVI_alpha
src_poststim.py     ← Post-stimulus Hilbert phase (poststim_ref_t) + ITC
src_lep.py          ← LEP feature extraction (N2/P2 amplitude, latency, N2-P2)
src_fooof.py        ← FOOOF / specparam aperiodic decomposition
src_write.py        ← CSV writers (trial-level, GA, FOOOF)
src_plot.py         ← 2-D QC figures (PSD, FOOOF, ERD, LEP, phase polar)
src_render.py       ← Offscreen 3-D brain renders (sLORETA maps, ROI scalar maps)

---

## 4. Module Reference

### Entry Points
| File | Role |
|---|---|
| `src_io.py` | `src_open_log`, `src_logmsg`, `src_close_log` - timestamped dual stdout/file logging. `src_find_set` - resolves the latest `_base.set` file for a subject with optional fallback search. `src_read_epochs` - MNE version-compatible EEGLAB epoch loader |
| `src_assets.py` | `src_load_fsaverage_assets` - validate the three required BEM files exist. `src_load_labels` - loads all parcellation labels for fsaverage. `src_build_custom_rois` - merges individual labels into named macro-ROIs; handles both `label-hemi` and `hemi-label` naming conventions and `BiHemiLabel` objects |
| `src_inverse.py` | `src_make_inverse_operator` — builds or reloads (cached) the forward solution, computes empirical noise covariance from the baseline window, and assembles the inverse operator. `src_apply_inverse_epochs` — applies sLORETA to every epoch, extracts ROI time courses; returns `(n_epochs, n_rois, n_times)` array. `src_apply_inverse_evoked` — applies sLORETA to the grand-average evoked response. |
| `src_spectral.py` | `src_psd_welch` — Welch PSD with adaptive window sizing for short epochs. `src_bandpower` — trapezoidal band integration. `src_cog` — spectral centre of gravity (PAF proxy). `src_bandpass_filter` — zero-phase Hamming FIR for Hilbert phase extraction. |
| `src_alpha_features.py` | `src_compute_alpha_features` — the 10 alpha metrics + psi_cog from a single PSD (mirrors `spec_compute_alpha_features_from_psd.m`). `src_compute_window_alpha_features` / `src_compute_ga_window_alpha_features` — same, looped over every (trial, ROI) for one window (whole / pre / post); also returns `psd_by_roi` (GA) for reuse by FOOOF. |
| `src_erd.py` | `src_compute_noise_stats` — per-ROI noise floor (eps0) + 5th-percentile pre-stim power thresholds, pooled across trials within each ROI (see module docstring for why this is per-ROI rather than pooled across ROIs, unlike MATLAB's pooling across channels). `src_compute_erd_metrics` — erd_slow, erd_fast, erd_pow_alpha_total, delta_erd, p5_flag; always numeric (epsilon floor), never NaN. |
| `src_compute_metric_deltas.py` | `src_compute_metric_deltas` — generic `delta_<name> = post_<name> - pre_<name>` for every field common to matched pre/post rows. |
| `src_add_prefix.py` | `src_add_prefix` / `src_add_prefix_rows` — prefix dict (or list-of-dict) keys with `whole_`/`pre_`/`post_`, skipping join columns. |
| `src_merge_rows.py` | `src_merge_rows` — merge any number of row-lists sharing key columns into one row per key; replaces the old pandas-merge-with-suffix-dedup pattern in `src_write.py`. |
| `src_prestim.py` | `src_compute_prestim_phase` / `src_compute_ga_prestim_phase` — Hilbert instantaneous phase at t=0. `src_compute_tvi_alpha` — normalised MSSD of the per-trial `pre_sf_balance` sequence (one scalar per subject × ROI); pre-stim-only by design. |
| `src_poststim.py` | `src_compute_poststim_phase` / `src_compute_ga_poststim_phase` — Hilbert instantaneous phase at `poststim_ref_t`. `src_compute_itc` — inter-trial phase coherence (ITC / PLV) over the post-stimulus window. |
| `src_lep.py` | `src_compute_lep_trial` — per-trial N2/P2 peak amplitude, latency, N2-P2 peak-to-peak, and window mean amplitudes. `src_compute_lep_ga` — same metrics applied to the grand-average time course. |
| `src_fooof.py` | `src_fit_fooof` — fits one PSD with specparam / fooof (compatibility shim handles both package names). Extracts aperiodic offset, exponent, optional knee, and the strongest alpha-band peak (CF, PW, BW). `src_compute_fooof_ga` — batch fit over all ROIs. |
| `src_write.py` | `src_write_trial_csv` / `src_write_ga_csv` — write an already-merged (via `src_merge_rows`) row list to CSV. **[V3.0.0]** No longer does its own pandas merging internally. `src_write_fooof_csv` — writes FOOOF GA metrics. |
| `src_plot.py` | `src_plot_ga_psd` — GA PSD overlay with alpha-band shading. `src_plot_fooof` — FOOOF model fit plot (works with both fooof and specparam APIs). `src_plot_erd` — per-trial erd_slow and erd_fast grid across ROIs. `src_plot_lep_ga` — GA LEP waveform per ROI with N2/P2 windows shaded. `src_plot_phase_polar` — polar histogram with mean resultant vector for pre- or post-stim phase. `src_plot_brain` — 2-D lateral brain snapshot using MNE's matplotlib backend. |
| `src_render.py` | `src_render_stc_timepoints` — offscreen 3-D renders of the GA sLORETA map at pre-stim mean, N2 peak, and P2 peak latencies. `src_render_roi_scalar` — paints per-ROI GA scalar metrics (e.g. `pre_sf_balance`, `delta_erd`) onto the fsaverage surface; signed metrics use a diverging RdBu_r scale, non-negative metrics use YlOrRd. |

---

## 5. Configuration Reference

All parameters live under `cfg.source` in the experiment JSON. Every key has a safe default; only three starred entries and required.

```json
{
    "exp": {
        "id": "exp01",
        "out_prefix": "26BB_62_",   // * REQUIRED - .set filename prefix
        "subjects": [1, 2, 3],      // * REQUIRED - subject IDs
        "out_dirname": "26ByBiosemi"  // recommeded; falls back to exp.id
    },

    "paths": {
        "da-root": "cifs/seminowicz/eegPainDatasets/CNED/da-analysis"
    },

    "source": {
        "enabled": true,

        "fsaverage": {
            "subjects_dir": "/path/to/mne/datasets"     // * REQUIRED - must contain fsaverage/
        },

        "input": {
            "stage_dir": "08_base",     // preprocessed .set stage to read from
            "allow_fallback_search": true   // retry without out_prefix if not found
        },

        "outputs": {
            "root": "AUTO"  // resolved to <da_root>/<exp_out>/source
        },

        "forward": {
            "mindist_mm": 5.0   // minimum source-to-skull distance
        },

        "inverse": {
            "method": "sLORETA",
            "snr": 3.0,     // signal-to-noise; lambda2 = 1/snr^2
            "loose": 0.2,   // orientation constraint (0=fixed, 1=free)
            "depth": 0.8,   // depth weighting exponent
            "pick_ori": null    // null=magnitude, "normal"=surface-normal
        },

        "noise_cov": {
            "tmin": -0.5,   // pre-stimulus analysis window (s)
            "tmax": 0.0
        },

        "poststim": {
            "tmin": 0.1,    // post-stimulus window; 0.1s avoids stimulus artefact
            "tmax": 0.8,
            "phase_ref_t": 0.2  // time at which post-stim Hilbert phase is sampled
        },

        "lep": {
            "n2_window": [0.15, 0.35],  // N2 peak search window (s)
            "p2_window": [0.25, 0.50]   // P2 peak search window (s)
        },

        "roi": {
            "parcellation": "aparc",    // Desikan-Killiany; or "aparc.a2009s"
            "mode": "mean_flip",    // ROI extraction mode
            "use_custom_rois": true,    // false = use all parcels
            "custom_rois": {
                "S1": ["postcentral-lh", "postcentral-rh"],
                "M1": ["precentral-lh", "precentral-rh"],
                "ACC": ["caudalanteriorcingulate-lh", "caudalanteriorcingulate-rh, "rostralanteriorcingulate-lh", "rostalanteriorcingulate-rh"],
                "Ins": ["insula-lh", "insula-rh"],
                "SII": ["supramarginal-lh", "supramarginal-rh"],
                "dlPFC": ["rostralmiddlefrontal-lh", "rostralmiddlefrontal-rh"]
            }
        },

        "spectral": {
            "fmin": 1.0,
            "fmax": 40.0,
            "alpha_band": [8.0, 12.0],
            "slow_alpha_band": [8.0, 10.0],
            "fast_alpha_band": [10.0, 12.0]
        },

        "fooof": {
            "enabled": true,
            "aperiodic_mode": "fixed",  // "fixed" or "knee"
            "peak_width_limits": [1.0, 12.0],
            "max_n_peaks": 6,
            "min_peak_height": 0.1,
            "peak_threshold": 2.0,
            "freq_range": [1.0, 40.0]
        },

        "qc": {
            "save_plots": true
        },

        "render": {
            "enabled": false,   // opt-in; requires pyvista
            "use_mesa": false,  // true on headless HPC without GPU
            "stc_enabled": true,    // render GA STC at pre/N2/P2
            "stc_clim_pct": [50, 99],
            "roi_enabled": true,    // render ROI scalar maps
            "roi_metrics": ["pre_sf_balance", "pre_paf_cog_hz", "delta_erd", "n2p2_amp"]
        }
    }
}
```

---

## 6. Metrics Extracted

### Alpha feature family — whole / pre / post / delta (per trial x ROI)

**[V3.0.0]** Every metric below is now computed for all three windows and
prefixed accordingly, plus a generic delta. Field names match
`spectral_core.m` exactly (same formulas, same eps0 = 1e-12) so channel-space
and source-space metrics can be compared directly by name.

| Metric family | `whole_<name>` (full epoch) | `pre_<name>` | `post_<name>` | `delta_<name>` = post − pre |
|---|---|---|---|---|
| Slow alpha [8–10 Hz] bandpower | `whole_pow_slow_alpha` | `pre_pow_slow_alpha` | `post_pow_slow_alpha` | `delta_pow_slow_alpha` |
| Fast alpha [10–12 Hz] bandpower | `whole_pow_fast_alpha` | `pre_pow_fast_alpha` | `post_pow_fast_alpha` | `delta_pow_fast_alpha` |
| Total alpha [8–12 Hz] bandpower | `whole_pow_alpha_total` | `pre_pow_alpha_total` | `post_pow_alpha_total` | `delta_pow_alpha_total` |
| Spectral centre of gravity (PAF proxy, Hz) | `whole_paf_cog_hz` | `pre_paf_cog_hz` (was `CoG_pre`) | `post_paf_cog_hz` | `delta_paf_cog_hz` |
| Relative slow alpha: slow / (alpha + ε) | `whole_rel_slow_alpha` | `pre_rel_slow_alpha` | `post_rel_slow_alpha` | `delta_rel_slow_alpha` |
| Relative fast alpha: fast / (alpha + ε) | `whole_rel_fast_alpha` | `pre_rel_fast_alpha` | `post_rel_fast_alpha` | `delta_rel_fast_alpha` |
| Slow/fast ratio: slow / (fast + ε) | `whole_sf_ratio` | `pre_sf_ratio` | `post_sf_ratio` | `delta_sf_ratio` |
| Log ratio: ln(slow+ε) − ln(fast+ε) | `whole_sf_logratio` | `pre_sf_logratio` (was `LR_pre`) | `post_sf_logratio` | `delta_sf_logratio` |
| Balance index: (slow−fast)/(slow+fast+ε) ∈ [−1,+1] | `whole_sf_balance` | `pre_sf_balance` (was `BI_pre`) | `post_sf_balance` | `delta_sf_balance` |
| Slow alpha fraction: slow/(slow+fast+ε) | `whole_slow_alpha_frac` | `pre_slow_alpha_frac` | `post_slow_alpha_frac` | `delta_slow_alpha_frac` |
| Interaction term: sf_balance × (paf_cog_hz − 10) | `whole_psi_cog` | `pre_psi_cog` (was `psi_cog`) | `post_psi_cog` | `delta_psi_cog` |

### ERD family (per trial x ROI)

Fractional pre→post change for the 3 power metrics only — not produced
generically above because a fractional change is not meaningful for a
bounded/ratio metric (sf_balance, sf_logratio, psi_cog, ...). Always a
numeric value using an epsilon floor (never NaN due to low pre-stim power);
low-power trials are flagged separately via `p5_flag` rather than dropped.

| Column | Description |
|---|---|
| `erd_slow` | Fractional slow-alpha desynchronization: (post − pre) / (pre + ε) |
| `erd_fast` | Fractional fast-alpha desynchronization: (post − pre) / (pre + ε) |
| `erd_pow_alpha_total` | Fractional total-alpha desynchronization: (post − pre) / (pre + ε) |
| `delta_erd` | `erd_slow − erd_fast` (< 0: fast-ERD dominant; > 0: slow-ERD dominant) |
| `p5_flag` | 1 if pre-stim slow or fast power is below the session's 5th percentile *for that ROI* (see `src_erd.py` for why this is per-ROI, not pooled across ROIs the way MATLAB pools across channels) |

### Phase (per trial x ROI)

Unchanged in scope — not part of the whole/pre/post/delta generalization.

| Column | Description |
|---|---|
| `slow_phase` | Hilbert instantaneous phase of 8–10 Hz signal at t = 0 (stimulus onset, radians) |
| `sin_phase` / `cos_phase` | Circular regressors for GAMMs |
| `slow_phase_post` | Hilbert phase of 8–10 Hz signal at `poststim_ref_t` (radians) |
| `sin_phase_post` / `cos_phase_post` | Circular regressors for GAMMs |

### Per subject x ROI (GA CSV only)

| Column | Description |
|---|---|
| `TVI_alpha` | Temporal variability index: normalized MSSD of the per-trial `pre_sf_balance` sequence. Pre-stim-only by design. |
| `itc_mean` | Mean inter-trial phase coherence over post-stim window |
| `itc_peak` | Peak ITC value |
| `itc_peak_latency_ms` | Latency of peak ITC (ms) |

GA versions of the whole/pre/post/delta/ERD families above also exist,
computed on the trial-averaged time course (see `src_alpha_features.py` /
`src_erd.py` docstrings for the exact GA convention, which mirrors
`spectral_core.m`'s split between PSD-averaged-first for `whole_*` and
feature-averaged-after for `pre_*`/`post_*`/`delta_*`). GA `erd_*` uses only
an epsilon guard, no `p5_flag` (the percentile of a single GA value is
degenerate).
 
### LEP features (per trial × ROI)
 
| Column | Description |
|---|---|
| `n2_amp` | N2 peak amplitude (most negative value in N2 window) |
| `n2_lat_ms` | N2 peak latency (ms) |
| `p2_amp` | P2 peak amplitude (most positive value in P2 window) |
| `p2_lat_ms` | P2 peak latency (ms) |
| `n2p2_amp` | N2-P2 peak-to-peak amplitude (P2 − N2) |
| `n2_mean` | Mean amplitude in N2 window |
| `p2_mean` | Mean amplitude in P2 window |
 
### FOOOF / aperiodic (per subject × ROI, FOOOF GA CSV)

**[V3.0.0]** Now fits the whole-epoch GA PSD (`whole_` window), not the
pre-stim GA PSD as in earlier versions — matches `spec_run_fooof_python.m`
/ MATLAB's FOOOF input.

| Column | Description |
|---|---|
| `fooof_offset` | Aperiodic offset (log-power intercept) |
| `fooof_exponent` | Aperiodic exponent (1/f slope) |
| `fooof_knee` | Knee frequency (NaN if `aperiodic_mode = "fixed"`) |
| `fooof_alpha_cf` | Alpha peak centre frequency (Hz); NaN if no peak detected |
| `fooof_alpha_pw` | Alpha peak power |
| `fooof_alpha_bw` | Alpha peak bandwidth (Hz) |

---

## 7. Output Layout

```
<da_root>/<exp_out>/source/
└── sub-XXX/
    ├── csv/
    │   ├── sub-XXX_source_trial.csv        per-trial metrics (all metrics above)
    │   ├── sub-XXX_source_ga.csv           GA pre+post+LEP+TVI+ITC per ROI
    │   └── sub-XXX_source_ga_fooof.csv     FOOOF aperiodic per ROI
    ├── figures/
    │   ├── sub-XXX_source_GA_roi_psd.png
    │   ├── sub-XXX_source_GA_<ROI>_fooof.png
    │   ├── sub-XXX_source_trial_erd.png
    │   ├── sub-XXX_source_GA_lep.png
    │   ├── sub-XXX_source_phase_prestim.png
    │   ├── sub-XXX_source_phase_poststim.png
    │   └── renders/                        (if render.enabled = true)
    │       ├── sub-XXX_render_stc_prestim_mean.png
    │       ├── sub-XXX_render_stc_n2_peak.png
    │       ├── sub-XXX_render_stc_p2_peak.png
    │       └── sub-XXX_render_roi_<metric>.png
    ├── fwd/
    │   └── sub-XXX_fwd.fif                 cached forward solution
    └── logs/
        └── sub-XXX_source_<timestamp>.log
```

---

## 8. Dependencies

```
mne >= 1.0
numpy
pandas
scipy
matplotlib
specparam (preferred) or fooof (legacy)
pyvista (optional - only required if render.enabled = true)
pyvistaqt (optional - may be required on some platforms)
```

**fsaverage template:** Run once to download:
```python
import mne
fs_dir = mne.datasets.fetch_fsaverage(verbose=True)
print(fs_dir) # pass this path as cfg.source.fsaverage.subjects_dir
```

**Headless rendering (HPC):** If no GPU is available, set `render.use_mesa = true` in the JSON. This routes PyVista through Mesa's OSMesa software rasterizer via `PYOPENGL_PLATFORM=osmesa`. Slower, but requres no X server or GPU.

---

## 9. Run Sequence

```bash
# A single experiment, all subjects in JSON
python exp01_source.py

# Override subjects from CLI
python exp01_souce.py --subjects 1 2 3 10

# Point to a non-default config
pythin exp01_source.py --cfg /path/to/config/exp01.json

# Multiple experiments
python exp01_source.py
python exp02_source.py
# ...
```

The forward solution for wach subject is cached to `fwd/sub-XXX_fwd.fif` on first run and reloaded on subsequent runs, making re-runs after parameter changes much faster.

---

## 10. Adding a New Experiment

Copy `exp01_source.py` and change two lines:

```python
exp_id = "exp02"    # change this
cfg_path = ".../config/exp02.json"  # and this (or pass --cfg)
```

Everything else - the full pipeline, all helpers, all output structyre - is inherited automatically through `souce_default` and `source_core`.