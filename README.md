# **Pain EEG Preprocessing & Analysis Toolkit**
*A modular, reproducible MATLAB pipeline for preprocessing, source localization, and analysis of pain-related EEG data.*

---

## **Overview**

This repository provides a complete, configurable pipeline for:

1. **EEG Preprocessing**
   - 9 experiment-specific configurations
   - Unified, reproducible core architecture
2. **Source Localization**
   - Support for 32-, 62-, and 64-channel caps
   - sLORETA and DeepSIF workflows (WIP)
3. **Analysis Framework**
   - Alpha-band feature extraction
   - Slow/fast alpha modelling
   - Generalized Mixed Additive Models (GAMMs)
   - Pain-brain relationships in source space

This toolkit is designed for **scalable, transparent, and publishable neuroscience workflows**, with complete documentation and reproducibility standards.

---

## **Pipeline Architecture**
```
Raw EEG
  -> Preprocessing (per experiment)
    -> Cleaned Epochs
      -> Source Localization (per cap type)
        -> ROI time series
          -> Analysis (shared model)
            -> Statistics & figures
```

---

## **Repository Structure**
```
pain-eeg-pipeline/
│
├── README.md
├── LICENSE
├── matlab/
│   ├── preproc/
│   │   ├── preproc_core.m
│   │   ├── preproc_default.m
│   │   ├── exp01_preproc.m
│   │   ├── exp02_preproc.m
│   │   ├── ...
│   │   ├── exp09_preproc.m
│   │   └── README.md
│   │
│   ├── source/
│   │   ├── source_core.m
│   │   ├── cfg_cap_20.m
│   │   ├── cfg_cap_32.m
│   │   ├── cfg_cap_64.m
│   │   ├── run_source_exp01.m
│   │   ├── run_source_exp02.m
│   │   └── README.md
│   │
│   ├── analysis/
│   │   ├── run_analysis.m
│   │   ├── alpha_models.m
│   │   ├── gamm_models.m
│   │   ├── stats_utils.m
│   │   └── README.md
│   │
│   ├── utils/
│   │   ├── load_cfg.m
│   │   ├── eegplot_wrapper.m
│   │   ├── file_utils.m
│   │   └── ...
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
├── examples/
│   ├── example_preproc_output/
│   ├── example_source_output/
│   ├── example_analysis_results/
│   └── notebooks/
│       └── validation.ipynb
│
└── config/
    ├── paths_template.json
    ├── exp01.json
    ├── exp02.json
    └── ...
```

---

## **Usage**

### **1. Preprocessing**
Each experiment has a wrapper script:
```matlab
matlab/preproc/exp01_preproc.m
```

Run:
```matlab
exp01_preproc
```

---

### **2. Source Localization**

```matlab
cfg = cfg_cap_32();
source_core(cfg);
```

Outputs standardized source-space signals.

---

### **3. Analysis**

After preprocessing and source localization:

```matlab
run_analysis
```

This computes:
- alpha-band metrics
- slow/fast alpha ratios
- GAMMs
- Statistical Outputs

---

## **Documentation**

See `/docs/` for detailed specifications:

- **Pipeline_Overview.md** - conceptual summary
- **Data_Standards.md** - naming conventions & directory schemas
- **Preprocessing_Design.md** - filters, ICA, referencing, epoching
- **Source_Localization_Design.md** - cap models, DeepSIF, sLORETA
- **Analysis_Methods.md** - alpha composition, GAMMs, stats

The **GtHub Wiki** provides tutorials and high-level explanations.

---

## **For Researchers**

This toolkit supports:
- EEG preprocessing using modern standards
- Source localization for pain-related paradigms
- Nonlinear statistical modelling
- High-reproducibility research workflows

It is intended for both new learners and experienced researchers.

---

## **Contact**

For questions or collaborations:
**darsenea@uwo.ca || drydena18.github.io || https://www.linkedin.com/in/dryden-arseneau/**
