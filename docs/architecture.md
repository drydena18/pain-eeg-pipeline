# Architecture

## 1. Layered Structure
This preprocessing framework is intentionally layered.  

Each file has one responsibility:  

| Layer | File | Responsibility |
|-------|------|----------------|
| Configuration | expXX.json | Define experiment-specific behaviour |
| Loader | load_cfg.m | Parse JSON |
| Paths | config_paths.m | Resolve I/O locations |
| Entrypoint | expXX_preproc.m | Select experiment + subjects |
| Normalizer | preproc_default.m | Validate + complete config |
| Engine | preproc_core.m | Execute preprocessing |
| Utilities | helpers | QC, metrics, prompting |

No file performs more than one conceptual task.

---

## 2. Experiment Registry

Experiments are declared centrally in:  
```
experiment_registry()
```
Each experiment specifies:
- Raw data folder
- Output folder name
- Output filename prefix

The pipeline engine does not need to be modified for new experiments.

---

# 3. Config-Driven Behaviour

All preprocessing decisions are defined in JSON.  

No filtering parameters are hard-coded in MATLAB.  

This allows:
- Version control of parameters
- Transparent documentation of decisions
- Easy experiment replication

---

## 4. Stage Execution Model
Each stage:
1. Checks whether output already exists
2. If yes → load and continue
3. If no → compute and save

This enables:
- Crash recovery
- Partial reruns
- Debugging individual steps
- Safe experimentation

---

## 5. ICA Training Design
ICA is trained on a temporary copy of the dataset to improve stability.  
Possible modifications to the training copy:
- Segment removal
- High-pass filtering
- Rank adjustment

After training, weights are transferred to the broadband dataset.  

This ensures:
- Stable decomposition
- Clean application to full data
- Clear seperation between training and analysis data