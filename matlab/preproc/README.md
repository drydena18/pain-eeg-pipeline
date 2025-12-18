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
