# Filesystem Layout
All outputs live under:
```
PROJ_ROOT/preproc/
```

---

## Per-Subject Structure
```
sub-XXX/  
01_filter/  
02_notch/  
03_resample/  
04_reref/  
05_initrej/  
06_ica/  
07_epoch/  
08_base/  
LOGS/  
QC/ 
```

---

## Directory Roles
| Folder | Purpose |
|--------|---------|
| 01_filter | High/low-pass outputs |
| 02_notch | Notch-filtered outputs |
| 03_resample | Resampled data |
| 04_reref | Re-referenced data |
| 05_initrej | Channel QC stage |
| 06_ica | ICA-decomposed data |
| 07_epoch | Epoched data |
| 08_base | Baseline corrected data |
| LOGS | Logs + QC figures |
| QC | Structued IC QC packets |

---

## Design Principles

- Raw data is read-only
- No stage overwrites previous outputs
- Each subject is fully isolated
- QC artifacts are seperated from EEG stage data
- Logging is centralied per subject

This layout supports:
- Reproducibility
- Resume logic
- Parallel execution
- Clear audit trails