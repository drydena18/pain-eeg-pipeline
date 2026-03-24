# Preprocessing Stages
This page documents every stage executed by <code>preproc_core.m</code>, its inputs, outputs, config parameters, default values, and resume behaviour.

---

## Overview
Stages run in a fixed order. Each stage is independently enabled/disabled in the JSON config. Completed stages are detected by the presence of their output <code>.set</code> file and skipped on rerun.

| # | Stage | Dir | Default Tag |
|---|-------|-----|-------------|
– | Import + Montage | *(no save)* | – |
1 | Filter | <code>01_filter</code> | <code>fir</code> |
2 | Notch | <code>02_notch</code> | <code>notch60</code> |
3 | Resample | <code>03_resample</code> | <code>rs500</code> |
4 | Re-Reference | <code>04_reref</code> | <code>reref</code> |
5 | INITREJ | <code>05_initrej</code> | <code>initrej</code> |
6 | ICA (+ICLabel) | <code>06_ica</code> | <code>ica</code> |
7 | Epoch | <code>07_epoch</code> | <code>epoch</code> |
8 | Baseline Correction | <code>08_base</code> | <code>base</code> |

---

## Pre-Stage: Import and Montage

Not saved. Runs before the first stage on every execution.

Raw file resolution (<code>resolve_raw_file.m</code>)  
Three resolution strategies are tried in order:
1. **BIDS Path** – <code>P.INPUT.EXP/sub-XXX/eeg/sub-XXX_task-<task>_eeg.bdf</code>
2. **JSON Pattern** – <code>cfg.exp.raw.pattern</code> is <code>sprintf</code>-formatted with <code>subjid</code> and resolved relative to <code>P.INPUT.EXP</code>
3. **Recursive Search** – scans <code>P.INPUT.EXP/&ast;&ast;</code> for <code>&ast;sub-XXX&ast;.bdf/eeg</code> (if <code>cfg.exp.raw.search_recursive = true</code>)

### Import
```matlab
EEG = pop_biosig(rawPath);
EEG = normalize_chan_labels(EEG); % strip spaces/dashes from labels
```

### Montage (optional)
Triggered by <code>cfg.exp.montage.enabled = true</code>.
1. Reads a CSV with columns <code>raw_label</code>, <code>std_label</code>
2. Optionally filters to A1-A32 and B1-B32 channels only (<code>select_ab_only</code>, default <code>true</code>)
3. Relabels channels; validates no duplicates introduced
4. Optionally applies coordinate lookup from <code>P.CORE.ELP_FILE</code> (<code>do_lookup</code>, default <code>true</code>)
5. Writes audit TSV: <code>LOGS/sub-XXX_channelmap_applied.tsv</code>

### Config
```json
"montage": {
    "enabled": true,
    "csv": "montage.csv",
    "select_ab_only": true,
    "do_lookup": true
}
```

### Deterministic RNG
```matlab
rng(double(subjid), 'twister')
```
Set before any stage runs. Ensures ICA produces identical results across reruns and machines for same subject.

---

## Stage 1: Filter
**Function**: <code>pop_eegfiltnew</code> (FIR)  
Applies a bandpass filter. Both highpass and lowpass are applied in one call.

### Config
```json
"filter": {
    "enabled": true,
    "tag": "fir",
    "type": "fir",
    "highpass_hz": 0.5,
    "lowpass_hz": 40.0
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
<code>highpass_hz</code> | <code>0.5</code> | High-pass cutoff (Hz) |
<code>lowpass_hz</code> | <code>40.0</code> | Low-pass cutoff (Hz) |
<code>type</code> | <code>"fir"</code> | Filter type label (informational; EEGLAB uses FIR via <code>pop_eegfiltnew</code>) |

**Output**: <code>01_filter/&lt;prefix&gt;&lt;subjid&gt;_fir.set</code>

---

## Stage 2: Notch
**Function**: <code>pop_eegfiltnew</code> with <code>revfilt = 1</code> (bandstop)  
Removes line noise. Implemented as a bandstop (notch) filter centered at <code>freq_hz</code> with half_bandwidth <code>bw_hz</code>.

### Config:
```json
"notch": {
    "enabled": true,
    "tag": "notch60",
    "freq_hz": 60,
    "bw_hz": 2
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
<code>freq_hz</code> | <code>60</code> | Notch center frequency (Hz) |
<code>bw_hz</code> | <code>2</code> | Half-bandwidth; stop bands is <code>[freq_hz - bw_hz, freq_hz + bw_hz]</code>

**Output**: <code>02_notch/&lt;prefix&gt;&lt;subjid&gt;_fir_notch60.set</code>

---

## Stage 3: Resample
**Function**: <code>pop_resample</code>  
Resamples the data to a target sampling rate. If the data is already at <code>target_hz</code>, this stage is skipped (no file saved, no tag appended)  
### Config:
```json
"resample": {
    "enabled": false,
    "tag": "rs500",
    "target_hz": 500
}
```

| Parameter | Default | Description |
|-----------|---------|-------------|
<code>target_hz</code> | <code>[]</code> |  Target sampling rate (Hz). Must be set if enabled.

**Output**: <code>03_resample/&lt;prefix&gt;&lt;subjid&gt;_fir_notch60_rs500.set</code>

---

