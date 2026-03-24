# Architecture
This page describes the structural design of the CNED EEG preprocessing framework: how the layers relate to each other, what each layer owns, and how configuration flows through the system.

---

## Layered Call Chain
The framework is organized into strict layers. Each layer has a single responsibility and delegates downward.

```
┌─────────────────────────────────────────────────────────┐
│  expXX_preproc.m          Experiment entrypoint         │
│  · Adds EEGLAB to path                                  │
│  · Loads JSON config                                    │
│  · Calls config_paths() → preproc_default()             │
└───────────────────────┬─────────────────────────────────┘
                        │
┌───────────────────────▼─────────────────────────────────┐
│  preproc_default.m        Config normalization layer    │
│  · Validates required fields                            │
│  · Resolves subjects (override > JSON > TSV)            │
│  · Fills all defaults (normalize_preproc_defaults)      │
│  · Validates epoch event_types are non-empty            │
│  · Calls preproc_core()                                 │
└───────────────────────┬─────────────────────────────────┘
                        │
┌───────────────────────▼─────────────────────────────────┐
│  preproc_core.m           Execution engine              │
│  · Loops over subjects                                  │
│  · Creates per-subject folder tree                      │
│  · Opens per-subject log file                           │
│  · Executes enabled stages in order                     │
│  · Each stage: check cache → process → save             │
└─────────────────────────────────────────────────────────┘
```
The spectral pipeline mirrors this exactly:
```
expXX_spectral.m -> spectral_default.m -> spectral_core.m
```

---

## Layer Responsibilities
<code>expXX_preproc / expXX_spectral.m</code>
</br>

The experiment entrypoint. Responsibilities.
- Add EEGLAB (and optionally helpers) to the MATLAB path
- Load the JSON config via <code>load_cfg</code>
- Resolve paths via <code>config_paths(exp_id, cfg)</code>
- Pass everything down to the default/normalization layer
- Accept an optional <code>subjects_override</code> argument

This is the **only** file that is experiment-specific. All logic lives downstream.

<code>config_paths.m</code>
</br>

Centralized path resolution. Responsibilities:
- Look up experiments in <code>experiments_registry()</code>
- Construct <code>P.RAW_ROOT</code>, <code>P.PROJ_ROOT</code>, <code>P.INPUT.EXP</code>, <code>P.RUN_ROOT</code>, <code>P.RESOURCE</code>
- Expose <code>P.CORE.ELP_FILE</code> and <code>P.CORE.PARTICIPANTS_TSV</code>
- Expose <code>P.NAMING.fname()</code> – the single naming function for all .set files
- Allow JSON overrides for <code>paths.raw_root</code> and <code>paths.proj_root</code>
- Create output directories (never RAW_ROOT)
</br>
</br>

<code>experiment_registry()</code>
A simple lookup table for all experiments. Each entry defines:  
| Field | Purpose|
|-------|--------|
<code>id</code> | Experiment identifier string
<code>raw_dirname</code> | Subfolder under <code>RAW_ROOT</code> containing <code>.bdf</code>/<code>.eeg</code> files
<code>out_dirname</code> | Subfolder under <code>PROJ_ROOT</code> for all outputs
<code>out_prefix</code> | Filename prefix for all <code>.set</code> outputs (e.g., <code>26BB_62_</code>)

To add a new experiment, add one entry here. No other files needs to change.
</br>
</br>

<code>preproc_default.m</code> / <code>spectral_default.m</code>
</br>

Config normalization and validation. Responsibilities:
- Assert required top-level fields exist (<code>cfg.exp</code>, <code>cfg.preproc</code>)
- Resolve the subject list (override argument -> <code>cfg.exp.subjects</code> -> <code>participants.tsv</code>)
- Call <code>normalize_preproc_defaults()</code> to fill all missing preproc fields
- Validate constraints (e.g., <code>epoch.event_types</code> must be non-empty if epoch is enabled)
- Print a run header to stdout
- Dispatch to the core engine
</br>
</br>

<code>preproc_core.m</code> / <code>spectral_code.m</code>
</br>

The execution engine. Responsibilities:
- Iterate over subjects
- Create the per-subject folder tree
- Open the per-subject log file
- For each enabled stage: check whether a cached output already exists (<code>maybe_load_stage</code>), skip if so, otherwise execute and save (<code>save_stage</code>)
- Carry cumulative <code>tags{}</code> forward – each completed stage appends its tag, which forms the output filename

---

## Configuration Flow
```
expXX.json
    │
    ├── load_cfg()              → raw struct from jsondecode()
    │
    ├── config_paths()          → P struct (paths, naming, registry entry)
    │
    └── preproc_default()
            │
            ├── normalize_preproc_defaults()   fills cfg.preproc.*
            ├── normalize_subject_ids()         normalizes subs vector
            │
            └── preproc_core()
                    │
                    └── per-stage blocks read from cfg.preproc.<stage>.*
```
Every parameter that drives behaviour comes from <code>cfg</code>. The core engine never harcodes values.

---

## File Naming Convention
All <code>.set</code> files follow the pattern:
```
<out_prefix><subjid_padded>_<tag1>_<tag2>..._<tagN>.set
```
For example:
```
26BB_62_001_fir_notch60_reref_initrej_ica.set
```
Tags are appended cumulatively as stages complete. This means the filename of the final output encodes the full preprocessing history. The naming function is <code>P.NAMING.fname(subjid, tags, prefix)</code>, implemented in <code>local_fname.m</code>.

---

## Stage Resumption Pattern
Each stage follows this identical pattern in <code>preproc_core.m</code>.:
```matlab
nextTag = char(string(cfg.preproc.<stage>.<tag>));
[EEG, tags, didLoad] = maybe_load_stage(stageDir, P, subjid, tags, nextTag, logf, EEG);
if ~didLoad
    % ... do the work ...
    tags{end+1} = nextTag;
    save_stage(stageDir, P, subjid, tags, EEG, logf);
end
```
<code>maybe_load_stage</code> checks whether a <code>.set</code> file with the expected name already exists. If it does, it loads it and sets <code>didLoad = true</code>, skipping all processing. This means you can kill and rerun a pipeline at any point – it resumes from the last completed stage.

---

## Helper Function Organization
Helper functions exist in two forms:
</br>

**Standalone <code>.m</code> files** (one function per file, e.g., <code>ensure_dir.m</code>, <code>logmsg.m</code>); these are on the MATLAB path and callable from anywhere.

<code>preproc_helpers.m</code> (namespace anchor): an older pattern (outdated) where multiple subfunctions live in one file. MATLAB resolves these by function name. Updated code should use standalone files.

<code>spec_*.m</code> files are spectral_specific helpers that mirror the preproc helpers but are prefixed to avoid name collisions.

---

## Multi-Experiment Design
All experiments share the same exectution engine. Adding a new experiment requires:

1. One entry in <code>experiment_registry()</code>
2. One JSON config file in <code>config/</code>
3. One entrypoint file (<code>expXX_preproc.m</code>) – copy and change the <code>exp_id</code> string

No logic changes are needed. The framework is designed so that experiment-specific decisions live entirely in JSON.