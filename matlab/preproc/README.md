# EEG Preprocessing Pipeline
**Experiment**: pain alpha dynamics  
**Author**: Dryden Arseneau  
**Labs**: Schabrun Lab, Seminowicz Lab  
**Data**: CNED (citation: Zhao, X., Zhou, J., Zhang, L. et al. A comprehensive EEG dataset of laser-evoked potentials for pain research. Sci Data 12, 1536 (2025). https://doi.org/10.1038/s41597-025-05900-1)

---

# 1. Purpose & Design Philosophy

This preprocessing pipeline was built to solve a specific problem:

> **How do we preprocess large, heterogenous EEG datasets in a way that is reproducible, inspectable, resumable, and defensible under peer review?**

Rather than prioritizing brevity, this pipeline prioritizes:

- Explicit decisions
- Human-in-the-loop quality control
- Full audit trails
- Config-driven reproducibility
- Stage-wide resumption
- Clear separation of responsibilites

This is **not** a one-off script.
It is a **preprocessing system**.

---

## 2. High-Level Architecture

This pipeline is intentionally layered:

```bash
expXX_preproc.m      ← experiment entrypoint
↓
preproc_default.m    ← config normalization & validation
↓
preproc_core.m       ← execution engine
↓
helper functions     ← QC, logging, plotting, prompting
```

Each layer has **one job** and does it explicitly.

---

## 3. File Responsibilities

### 3.1 'expXX_preproc.m' - Entrypoint

**Purpose:**  
Defines **what experiment is being run* and *which subjects to process*.  

**Responsibilities:**
- Load experiment paths
- Load JSON configuration
- Optionally override subject list
- Dispatch pipeline

This is the **only file you should normally edit** when running the pipeline.

#### Typical Usage
```matlab
exp01_preproc()
```

#### Run a single subject for debugging
```matlab
subjects_override(12);
preproc_default(exp_id, P, cfg, subjects_override(;
```

### 3.2 'preproc_default.m' - Configuration Normalization

**Purpose:**  
Transforms a flexible JSON configuration into a **validated, complete MATLAB struct**.  

This file answers:
> "Is this configuration safe and complete to execute?"

**Responsibilities:**  
- Validate required fields
- File defaults for optional fields
- Normalize subject IDs
- Normalize preprocessing blocks
- Enforce logical constraints

Does **not** perform EEG preprocessing.

### 3.3 'preproc_core.m' - Execution Engine  

**Purpose:**  
Executes preprocesing *exactly as specified* by the configuration.  

**Responsibilities:**
- Iterate over subjects
- Resolve raw data files
- Apply preprocessing stages
- Save intermediate outputs
- Generate QC artifacts
- Prompt for manual decisions
- Log all actions

Nothing irreversible occurs without either a log entry or a human decision.

---

## 4. Stage-Based Processing  

Each preprocessing step is treated as an **independant stage**.  

### 4.1 Why stages?  

Stages allow:
- Resumption after crashes
- Skipping completed steps
- Inspection of immediate data
- Debugging without rerunning everything

Each stage produces a tagged .set file:
```bash
26BB_64_012_fir_notch60_rs500_reref.set
```

### 4.2 Stage order
  1. Import raw .eeg
  2. Filter (high-pass / low-pass)
  3. Notch filter
  4. Resample
  5. Re-reference
  6. INITReJ (channel / segment QC)
  7. ICA
  8. ICLabel
  9. Manual IC rejection
  10. Epoching
  11. Baseline correction

Each stage can be enabled or disabled via JSON.

---

## 5. Filesystem Layout

All outputs live under:
```bash
da-analysis/expXX/
```

**Stage Directories:**

| Folder      | Description              |
| ----------- | ------------------------
| _filter/    | High/low-pass filtering  |
| _resample/  | Resampled data           |
| _reref/     | Re-referenced data       |
| _initrej/   | Initial channel QC       |
| _ica/       | ICA-cleaned data         |
| _epoch/     | Epoched data             |
| _base/      | Baseline-corrected data  |
| logs/       | Logs and QC figures      |

---

## 6. Configuration via JSON

All preprocessing behaviour is controlled via expXX.json.

### Why JSON?
- Human-readable
- Version-controllable
- Experiment-specific
- Decoupled from MATLAB code

No preprocessing decitions are hard-coded.

---

## 7. Montage & Channel Handling

### 7.1 Biosemi Channel Mapping

BioSemi data uses A1-32 and B1-32 channel labels.

**This pipeline:**
- Selects only A/B channels
- Maps to standard 10-20 labels via CSV
- Checks for duplicate labels
- Warns about missing midline electrodes
- Optionally applies .elp coordinate lookup

### 7.1.1 Channel Mapping Audit Trail

Every montage application produces:
```bash
logs/sub-XXX_channelmap_applied.tsv
```

### 7.2 AntNeuro and BrainProducts to come...

---

## 8. INITREJ - Initial Channel Quality Control

INITREJ is intentionally **conservative and manual**.

### 8.1 Automated suggestions (non-destrictive)

Suggestions are based on:
- EEGLAB probability z-scores
- EEGLAB kurtosis z-scores
- Channel RMS and STD
- PSD-based metrics
    - line-noise ratio
    - high frequency ratio
    - drift ratio
    - alpha dominance

No channel is removed automatically.

### 8.2 QC Outputs

Saved to logs/:
- RMS / STD histograms
- RMS / STD bar plots
- Topoplots
- PSD overview (median + IQR)
- PSD overlays for suggested channels
- Channel metrics CSV

### 8.3 Manual Decision

The user is prompted to select channels to interpolate.  

Default: **none**  

All decisions are logged.

---

## 9. ICA & ICLabel

### 9.1 ICA Training

ICA can be trained on:
- full continuous data
- or a cleaned copy with large-amplitude segments removed

Segment removal is:
- threshold-based
- logged
- optional
- user-confirmed

### 9.2 ICLabel Integration

ICLabel is used **only for suggestions**.
- No auto-rejection
- Thresholds defined in JSON
- Suggested ICs are retained unless manually removed

### 9.3 IC QC Packets

For each suggested IC:
- scalp topography
- PSD
- activation time series

Saved to:
```bash
logs/sub-XXX_icqc/
```

### 9.4 Manual IC Rejection

The user manually chooses ICs to remove.  

Logged outcomes:
- suggested
- removed
- kept

---

## 10. Logging & Reproducibility

### 10.1 Per-Subject Logs

Each subject has:
```bash
logs/sub-XXX_preproc.log
```

Logs include:
- parameters
- warnings
- QC summaries
- manual decisions
- skipped/resumed stages

### 10.2 Deterministic RNG
```matlab
rng(subjid, 'twister')
```
Ensures reproducible ICA results.

---

## 11. Why This Pipeline Is Large

Many EEG pipelines are short because they:
- hide assumptions
- hard-code parameters
- overwrite data
- rely on GUI state
- silently fail

This pipeline is long because it:
- makes assumptions explicit
- forces human decisions
- never overwrites data
- logs everything
- scales to multi-experiment datasets

**Length here is a feature, not a flaw**.

---

## 12. Intended Audience

This pipeline is written for:
- the author (future-proofing)
- labmates
- reviewers
- committee members
- future trainees

If someone asks *"why did you do X?"* – the answer is in:
- the JSON
- the logs
- this README

---

## 13. Future Extensions

- QC enable/disable flags
- Headless batch mode
- Helper modulization
- Auto-generated methods summaries
