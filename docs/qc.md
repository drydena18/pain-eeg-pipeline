# Quality Control (QC)

## Philosophy

OC in this framework is:
- Conservative
- Human-in-the-loop
- Non-destructive

Automated metrics provide suggestions only.  
The user makes final decisions.

---

## INITREJ (Channel QC)

Automated metrics include:
- RMS
- Standard deviation
- Kurtosis
- Probability
- PSD metrics:
    - Line noise rario
    - High-frequency ratio
    - Drift ratio
    - Alpha dominance

Outputs:
- Histogram
- Topoplots
- PSD Overlays
- CSV Summaries

No channel is removed automatically.

---

## ICA QC
ICLabel is used for suggestions only.  

For each suggested IC:
- Scalp topography
- PSD
- Activation time series

Saved under:
```
sub-XXX/QC
```
Manual rejection is required.

---

## Logging
All QC decisions are recorded in:
```
sub-XXX/LOGS/sub-XXX_preproc.log
```
Suggested, removed, and retained components are explicitly documented.