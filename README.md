# **Pain EEG Preprocessing & Analysis Toolkit**
This reporsitory exists to standardize pain EEG processing and analysis while enabling transparent, shareable, and extensible methods for neuroscience research.
*A modular, reproducible pipeline for preprocessing, source localization, and analysis of pain-related EEG data.*

---

## **Overview**

This repository provides a complete, configurable pipeline for:

1. **EEG Preprocessing (MATLAB)**
   - 9 experiment-specific configurations
   - Unified, reproducible core architecture
2. **Source Localization (Python)**
   - Support for 32-, 62-, and 64-channel caps
   - sLORETA and DeepSIF workflows (WIP)
   - Standardized source-space outputs
3. **Analysis & Modelling (R & Python)**
   - Alpha-band feature extraction
   - Slow/fast alpha modelling
   - Generalized Mixed Additive Models (GAMMs)
   - Pain-brain relationships in source space

This toolkit is designed for **scalable, transparent, and publishable neuroscience workflows**, with complete documentation and reproducibility standards.

---

## **Key Features**

- **Multi-experiment preprocessing** with shared architecture
- **Automatic ICA cleaning workflows**
- **Cap-specific source localization pipelines**
- **DeepSIF-ready integration** (WIP)
- **Unified ROI-level outputs for all experiments**
- **Cross-language reproducibility** (MATLAB > Python > R)
- **Publication-ready figures + GAMM models**
- **Designed for open methods papers and co-authorship pipelines**

## **Pipeline Architecture**
```
Raw EEG
  -> Preprocessing (MATLAB; per experiment)
    -> Cleaned sensor-space data
      -> Source Localization (Python; per cap type)
        -> ROI-level source signals
          -> Analysis & statistics (R / Python; shared model)
            -> Figures, tables, reports
```

---

## **Repository Structure**
```
pain-eeg-pipeline/
│
├── README.md
├── LICENSE
│
├── matlab/                      # Preprocessing
│   ├── preproc/
│   │   ├── preproc_core.m
│   │   ├── preproc_default.m
│   │   ├── exp01_preproc.m
│   │   ├── exp02_preproc.m
│   │   ├── …
│   │   ├── exp09_preproc.m
│   │   └── README.md
│   │
│   ├── utils/
│   │   ├── load_cfg.m
│   │   ├── file_utils.m
│   │   └── …
│   └── README.md
│
├── python/                      # Source localization
│   ├── source/
│   │   ├── source_core.py
│   │   ├── cfg_cap_32.py
│   │   ├── cfg_cap_62.py
│   │   ├── cfg_cap_64.py
│   │   ├── run_source_exp01.py
│   │   ├── run_source_exp02.py
│   │   └── README.md
│   │
│   ├── utils/
│   │   ├── io_utils.py
│   │   ├── plotting_utils.py
│   │   └── …
│   └── requirements.txt
│
├── R/                           # Analysis & statistics
│   ├── analysis/
│   │   ├── run_analysis.R
│   │   ├── alpha_features.R
│   │   ├── gamm_models.R
│   │   ├── stats_utils.R
│   │   └── README.md
│   └── renv/ (optional)
│
├── config/
│   ├── paths_template.json
│   ├── exp01.json
│   ├── exp02.json
│   └── …
│
├── docs/
│   ├── Pipeline_Overview.md
│   ├── Data_Standards.md
│   ├── Experiment_Configs.md
│   ├── Preprocessing_Design.md
│   ├── Source_Localization_Design.md
│   ├── Analysis_Methods.md
│   ├── Figures/
│   └── Diagrams/
│
└── examples/
├── example_preproc_output/
├── example_source_output/
├── example_analysis_results/
└── notebooks/
└── validation.ipynb
```

---

## **Dependencies**

### MATLAB (Preprocessing)
- MATLAB R20XX+
- EEGLAB
- Signal Processing Toolbox (recommended)

### Python (Source Localization)
Installation:
```bash
cd python
pip install -r requirements.txt
```

### R (Analysis & GAMMs)
Recommended Packages:
- mgcv
- tidyverse
- data.table
- lme4 / glmmTMB
- renv (optional for reproducibility)

## **Usage**

### **1. Preprocessing**

```matlab
cd matlab/preproc
exp01_preproc
```
Outputs cleaned, epoched data in:
```bash
/derivatives/preproc/EXP01/sub-XX/
```


### **2. Source Localization**

```python
cd python/source
python run_source_exp01.py
```
Outputs standardized source-space data in:
```bash
/derivatives/source/EXP01/sub-XX/
```

### **3. Analysis**

```r
setwd('R/analysis')
source("run_analysis.R")
run_analysis()
```
Outputs statistical results & figures in:
```bash
/derivatives/analysis/
```
---

## **Versioning**
This toolkit follows semantic versioning:
- **v0.x** - Development
- **v1.x** - First stable full pipeline
- **v2.x** - DeepSIF integration + publication release

---

## **Documentation**

See `/docs/` for detailed specifications:

- **Pipeline_Overview.md** - conceptual summary
- **Data_Standards.md** - naming conventions & directory schemas
- **Preprocessing_Design.md** - filters, ICA, referencing, epoching
- **Source_Localization_Design.md** - cap models, DeepSIF, sLORETA
- **Analysis_Methods.md** - alpha composition, GAMMs, stats

The **GitHub Wiki** provides tutorials and high-level explanations.

---

## **For Researchers**

This toolkit supports:
- EEG preprocessing using modern standards
- Source localization for pain-related paradigms
- Nonlinear statistical modelling
- High-reproducibility research workflows

It is intended for both new learners and experienced researchers.

---

## **Contributing**

Contributions, issues, and feature requests are welcome.
Please open an issue or submit a pull request via GitHub.

---

## **Contact**

**Dryden Arseneau**
Email: darsenea@uwo.ca
Website: drydena18.github.io 
LinkedIn: https://www.linkedin.com/in/dryden-arseneau/**
