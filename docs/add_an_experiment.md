# Adding a New Experiment

This page walks through adding a new experiment to the pipeline from scratch. "New experiment" means a new dataset with its own raw data folder, its own output directory, and potentially its own channel layout or event structure.

Nothing in the core engine (`preproc_core.m`, `spectral_core.m`) needs to change. Everything experiment-specific lives in three places: the registry, a JSON config file, and an entrypoint `.m` file.

---

## Step 1: Register the Experiment

Open `experiment_registry.m` and add one entry to the returned struct:

```matlab
function R = experiment_registry()

R = struct();

% existing experiments
R.exp01 = struct( ...
    'id',          'exp01', ...
    'raw_dirname', '26ByBiosemi', ...
    'out_dirname', '26ByBiosemi', ...
    'out_prefix',  '26BB_62_' ...
);

% new experiment
R.exp02 = struct( ...
    'id',          'exp02', ...
    'raw_dirname', 'MyNewDataset', ...
    'out_dirname', 'MyNewDataset', ...
    'out_prefix',  'MND_' ...
);

end
```

| Field | What it controls |
|-------|-----------------|
| `id` | Key used to look up this experiment. Must be a valid MATLAB struct field name (no hyphens). |
| `raw_dirname` | Subfolder under `RAW_ROOT` that contains the raw `.bdf`/`.eeg` files |
| `out_dirname` | Subfolder under `PROJ_ROOT` for all outputs (preproc, spectral, resource) |
| `out_prefix` | Prefix for all `.set` output filenames (e.g. `MND_001_fir_notch60_...set`) |

The `id` field is also used as the task label in BIDS path construction if `cfg.exp.task` is not set.

---

## Step 2: Create the JSON Config

Create `config/exp02.json`. The only truly required fields are `exp.id`, `exp.out_prefix`, and `preproc.epoch.event_types` (if epoching is enabled).

Start with the minimal config and add only what differs from the defaults:

```json
{
    "exp": {
        "id":         "exp02",
        "out_prefix": "MND_"
    },
    "preproc": {
        "epoch": {
            "event_types": ["Stim"]
        }
    }
}
```

### Common additions

**Custom subject list** (instead of reading `participants.tsv`):
```json
"exp": {
    "subjects": [1, 2, 3, 10, 11]
}
```

**Custom raw filename pattern** (if not BIDS-organized):
```json
"exp": {
    "raw": {
        "pattern": "Subject_%03d.bdf"
    }
}
```

**Biosemi montage relabelling:**
```json
"exp": {
    "montage": {
        "enabled": true,
        "csv": "montage.csv",
        "select_ab_only": true,
        "do_lookup": true
    }
}
```

Place the montage CSV in `<out_dirname>/resource/montage.csv` (i.e. `PROJ_ROOT/MyNewDataset/resource/montage.csv`). The CSV needs two columns: `raw_label` and `std_label`.

**Different filter or reference settings:**
```json
"preproc": {
    "filter": {
        "highpass_hz": 1.0,
        "lowpass_hz":  45
    },
    "reref": {
        "mode": "channels",
        "channels": [65, 66]
    }
}
```

**Different epoch window:**
```json
"preproc": {
    "epoch": {
        "event_types": ["S1", "S2"],
        "tmin_sec": -2.0,
        "tmax_sec":  3.0
    },
    "baseline": {
        "window_sec": [-1.5, -0.5]
    }
}
```

**Running on a different machine** (different paths):
```json
"paths": {
    "raw_root":  "/Volumes/lab/CNED",
    "proj_root": "/Users/me/analysis"
}
```

See [Configuration Reference](configuration.md) for every available field.

---

## Step 3: Create the Entrypoint File

Copy `exp01_preproc.m` to `exp02_preproc.m` and change exactly one line:

```matlab
function exp02_preproc(subjects_override)

exp_id = "exp02";   % ← change this

% everything else stays identical
```

If EEGLAB lives in a different path on this machine, also update the `addpath` line.

---

## Step 4: Place the Resource Files

The pipeline looks for shared resources in `PROJ_ROOT/<out_dirname>/resource/`. Create that folder and place any needed files there.

The directory is created automatically when `config_paths` runs, but the files themselves must be placed manually.

| File | Required when |
|------|--------------|
| `standard-10-5-cap385.elp` | `montage.do_lookup = true` or `channel_locs.use_elp = true` |
| `montage.csv` | `montage.enabled = true` |
| `participants.tsv` | `cfg.exp.subjects` not specified in JSON |

`participants.tsv` must have at least one of these column names: `participant_id`, `subjid`, or `subject`. Subject IDs can be integers or strings like `sub-001`.

---

## Step 5: Run

```matlab
% Run all subjects from participants.tsv
exp02_preproc();

% Run specific subjects
exp02_preproc([1 2 3]);
```

On first run the pipeline creates the full output folder tree under `PROJ_ROOT/MyNewDataset/preproc/` and begins processing. If it is interrupted, rerun the same command — completed stages are detected automatically and skipped.

---

## Step 6: Add the Spectral Pipeline (Optional)

If you want spectral features for this experiment, repeat steps 3 and 4 for the spectral side.

Copy `exp01_spectral.m` to `exp02_spectral.m` and change `exp_id = "exp02"`.

Add a `"spectral"` block to `config/exp02.json`:

```json
"spectral": {
    "enabled": true,
    "input_stage": "08_base",

    "psd": {
        "fmin_hz": 1,
        "fmax_hz": 80,
        "window_sec": 2.0,
        "overlap_frac": 0.5
    },

    "alpha": {
        "alpha_hz": [8, 12],
        "slow_hz":  [8, 10],
        "fast_hz":  [10, 12]
    },

    "fooof": {
        "enabled": true,
        "python_exe": "python3",
        "script_path": "/path/to/fooof_bridge.py",
        "fmin_hz": 2,
        "fmax_hz": 40
    }
}
```

Run:
```matlab
exp02_spectral();
exp02_spectral([1 2 3]);
```

---

## Checklist

```
□ Added entry to experiment_registry.m
□ Created config/exp02.json
    □ exp.id matches registry key
    □ exp.out_prefix set
    □ preproc.epoch.event_types set (if epoching)
□ Created exp02_preproc.m (exp_id updated)
□ Created exp02_spectral.m (exp_id updated, if needed)
□ Placed resource files in PROJ_ROOT/<out_dirname>/resource/
    □ participants.tsv (if subjects not in JSON)
    □ montage.csv (if montage.enabled)
    □ standard-10-5-cap385.elp (if do_lookup or use_elp)
□ Verified raw data is accessible at RAW_ROOT/<raw_dirname>/
□ Test run: exp02_preproc([1])
```

---

## Troubleshooting

**`Unknown exp_id "exp02"` error**
The registry lookup strips underscores from the `exp_id` before searching. Make sure the struct field name in `experiment_registry.m` matches — `exp02` not `exp_02`.

**`Raw file not found for sub-001. Skipping.`**
The three-step file resolution (BIDS → pattern → recursive search) failed. Check that `raw_dirname` in the registry matches the actual folder name under `RAW_ROOT`, and that the subject folder or raw file exists. Add `cfg.exp.raw.pattern` if the filename doesn't follow BIDS naming.

**`cfg.exp.subjects empty and participants.tsv not found`**
Either add `"subjects": [1, 2, 3]` to the JSON, or place a `participants.tsv` in `PROJ_ROOT/<out_dirname>/resource/` or in `RAW_ROOT/<raw_dirname>/`.

**`Montage CSV not found`**
The CSV path is resolved relative to `P.RESOURCE` (`PROJ_ROOT/<out_dirname>/resource/`). Place `montage.csv` there, or use an absolute path in `cfg.exp.montage.csv`.

**`Montage CSV expects channel not present in EEG: A33`**
The CSV contains a `raw_label` that doesn't exist in the loaded EEG. This usually means `select_ab_only` is `true` (keeping only A1–A32 and B1–B32) but the CSV lists channels outside that range, or the channel labels didn't normalize correctly. Check the debug output in the log for the actual labels present before relabelling.

**`cfg.preproc.epoch.enabled is true, but epoch.event_types is empty`**
Add `"event_types"` to the `epoch` block in your JSON. Run `validate_events_before_epoch` mentally first: check what event types your raw data actually contains (look at `EEG.event` after import in an interactive session) before choosing.