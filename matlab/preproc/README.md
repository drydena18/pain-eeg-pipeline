# CNED EEG Preprocessing Framework (MATLAB)

**Author:** Dryden Arseneau  
**Affiliations:** Schabrun Lab · Seminowicz Lab  
**Dataset:** CNED – Zhao et al., Sci Data (2025)  
https://doi.org/10.1038/s41597-025-05900-1

---

# 1. Overview

This repository contains a **config-driven, stage-based EEG preprocessing framework** designed for large, heterogenous, multi-experiment pain EEG datasets.  
The framework was built to address the following problem:  
> How can we preprocess large EEG datasets in a way that is reproducible, inspectable, resumable, and defensible under peer review?

This is not a one-off script.  
It is a modular preprocessing system.

---

# 2. Design Principles

The framework prioritizes:
- Explicit preprocessing decisions
- Human-in-the-loop quality control
- Full audit logging
- JSON-driven reproducibility
- Stage-based resumption
- Strict separation of responsibilites
- Deterministic ICA training

Raw data is never modified.   
All outputs are isolated per subject and per stage.

---

# 3. High-Level Architecture
```
expXX.json          ← Experiment configuration  
↓  
load_cfg.m          ← JSON loader  
↓  
config_paths.m      ← Path resolution & experiment registry  
↓  
expXX_preproc.m     ← Experiment extrypoint  
↓  
preproc_default.m   ← Config normalization & validation  
↓  
preproc_core.m      ← Execution engine  
↓  
Helper functions    ← QC, logging, prompting, metrics  
```

Each layer has a single responsibility.

---

# 4. Multi-Experiment Support
Experiments are registered in:
```
config_paths.m → experiment_registry()
```

Each experiment defines:
- Raw data directory
- Output directory name
- Filename prefix

All experiments use the same execution engine (`preproc_core.m`).   
This ensures consistent preprocessing across datasets. 

---

# 5. Stage-Based Processing
Preprocessing is organized into independent stages:  
1. Filter
2. Notch
3. Resample
4. Re-reference
5. INITREJ (channel QC)
6. ICA
7. ICLabel
8. Epoch
9. Baseline  

Each stage:
- Saves a tagged `.set` file
- Writes to its own directory
- Can be resumed independently

---

# 6. Reproducibility

### Deterministic ICA
```MATLAB
rng(subjid, 'twister')
```
ICA results are reproducible per subject across machines and reruns.  

### Logging
Each subject generates:
```
sub-XXX/LOGS/sub-XXX_preproc.log
```
Logs include:
- Parameters
- QC summaries
- Manual decisions
- Warnings
- Resume events

---

# 7. Intended Audience
This framework is written for:
- The author (future-proofing)
- Lab members
- Reviewers
- Committee members
- Future trainees

---

# 8. Documentation
Detailed technical documentation is available in:  
docs/  
architecture.md  
stages.md  
qc.md  
filesystem.md  
