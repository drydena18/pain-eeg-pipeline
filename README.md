# CNED EEG Preprocessing & Analysis Framework (MATLAB/R/PYTHON)

**Author**: Dryden Arseneau  
**Affiliations**: Human Pain Discovery Lab · Seminowicz Pain Imaging Lab · Univeristy of Western Ontario, Neuroscience, Schulich School of Medicine and Dentistry  
**Dataset**: CNED – Zhao et al., _Sci Data_ (2025) – <https://doi.org/10.1038/s41597-025-05900-1>  </br>
[![DOI](https://zenodo.org/badge/19358703.svg)](https://doi.org/10.5281/zenodo.19358703)

---

## Overview

This repository contains a **config-driven, stage-based EEG preprocessing and analysis framework** designed for large, heterogeneous, multi-experiment pain EEG datasets. It was built to answer one question:
> How do we preprocess large EEG datasets in a way that in reproducible, inspectable, resumable, and defensible under peer review?  

The framework is not a one-off scrip. It is a modular system spanning raw EEG import through spectral feature extraction and statistical modelling.  

---

## Repository Structure

```
pain-eeg-pipeline/
├── matlab/
│   ├── preproc/
│   │   ├── expXX_preproc.m
│   │   ├── preproc_default.m
│   │   ├── preproc_core.m
│   ├── utils/
│   │   ├── config_paths.m/      Centralized path resolution + experiment registry
│   ├── helpers/        Individual standalone helper functions
│   ├── spectral/
│   │   ├── expXX_spectral.m
│   │   ├── spectral_default.m
│   │   ├── spectral_core.m
│   │   ├── helpers/      Individual standalone helper functions
│   │   │   ├── spec_*.m      All spectral helpers prefixed with spec_
├── 
├── python/
│   ├── spectral/
│   │   ├── fooof_bridge.py     CLI  bridge: MATLAB -> specparam (FOOOF)
├── 
├── R/
│   ├── analysis/
│   │   ├── behavioural-analysis/
│   │   │   ├── merge_behavioural.R   Merge raw behavioural CSVs -> master
│   │   │   ├── merge_participants_into_behavioural.R    Add demographics + cap size
│   │   │   ├── merge_spectral_behaviour.R     Merge EEG spectral features -> master
│   │   │   ├── run_gamm_alpha_metrics.R      First-pass GAMM workflow
│   │   │   ├── run_gamm_alpha_metrics_v2.R    v2 GAMM with aperiodic controls + tensors
│   │   ├── experiment/    Experiment specific trial-by-trial behavioural data
```

---

## Design Principles

- **Explicit decisions** – every parameter lives in a JSON config; nothing is hardcoded
- **Human-in-the-loop QC** – bad channel interpolation and IC rejection are manual prompts backed by automated suggestions
- **Full audit logging** – every run produces a timestamped .log per subject
- **JSON-driven reproducibility** – rerunning the same JSON and raw data produces identical results
- **Stage-based resumption** – each preprocessing stage saves a tagged .set file; completed stages are skipped on rerun
- **Strict separation of responsibilities** – each layer of the call chain has a single job
- **Deterministic ICA** – `rng(subjid, 'twister')` seeds the RNG per subject before ICA
- **Raw data is never modified** – all outputs are isolated per subject and per stage under `PROJ_ROOT`

---

## Call Chains

### Preprocessing

```
pain-eeg-pipeline/
├── expXX_preproc.m
│   ├── preproc_default.m (validate + normalize config, resolve subjects)
│   │   ├── preproc_core.m (execute stage loop per subject)
```

### Spectral

```
pain-eeg-pipeline/
├── expXX_spectral.m
│   ├── spectral_default.m (validate + normalize config, resolve subjects)
│   │   ├── spectral_core.m (trial-wise PSD,  features, FOOOF, CSV, plots)
```

### R Analysis

```
pain-eeg-pipeline/
├── merge_behavioural.R
│   ├── merge_participants_into_behaviour.R
│   │   ├── merge_spectral_behaviour.R
│   │   │   ├── run_gamm_alpha_metrics.R / run_gamm_alpha_metrics_v2.R
```

---

## Quick Start

### 1. Prerequisites

- MATLAB (R2021b or later recommended)
- EEGLAB with BIOSIG and ICLabel plugins
- Python3 with fooof / specparam installed (for FOOOF step only)
- R with mgcv, readr, dplyr, ggplot2, purrr

### 2. Add a new experiment

Register it in `config_paths.m/`

```matlab
R.exp02 = struct(
  'id', 'exp02', ...
  'raw_dirname', 'MyRawFolder', ...
  'out_dirname', 'MyOutputFolder', ...
  'out_prefix', 'EXP02_'
);
```

### 3. Configure

Create `utils/exp02.json`:

```json
{
  "exp": {
    "id": "exp02",
    "out_prefix": "EXP02_",
    "subjects": [1, 2, 3]
  },

  "preproc": {
    "filter": {
      "enabled": true,
      "highpass_hz": 0.5,
      "lowpass_hz": 40
    },
    "notch":{
      "enabled": true,
      "freq_hz": 60
    },
    "resample": {
      "enabled": false
    },
    "reref": {
      "enabled": true,
      "mode": "average"
    },
    "initrej": {
      "enabled": true
    },
    "ica": {
      "enabled": true,
      "method": "runica",
      "iclabel": {
        "enabled": true,
        "thresholds": {
          "eye": 0.8,
          "muscle": 0.8
        }
      }
    },
    "epoch": {
      "enabled": true,
      "event_types": ["S1"],
      "tmin_sec": -1.0,
      "tmax_sec": 2.0
    },
    "baseline": {
      "enabled": true,
      "window_sec": [-0.5, 0]
    }
  }
}
```

### 4. Run

```matlab
exp01_preproc(); % All subjects from participants.tsv
exp01_preproc([1 2 3]); % Override subject list
```

---

## Further Documentation

| Page | Contents |
|------|----------|
Architecture | Call chain, layer responsibilities, config flow
Filesystem | Full directory tree, file naming, what lives where
Stages | Every preprocessing stage; inputs, outputs, parameters, resume logic
QC | QC outputs, plots, logs, manual decision plots
Spectral | Trial-wise spectral pipeline, features, FOOOF bridge
R Analysis | Behavioural merging, demographic merging, GAMM workflow
