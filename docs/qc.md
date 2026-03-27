# Quality Control (QC)

This page documents all QC outputs, the information each one contains, where to find them, and the manual decision points in the pipeline.

---

## Philisophy
The framework uses a **suggest-then-confirm** model for all QC decisions. Automated methods flag candidates, a human operator makes the final call. Every decision, including the choice to accept a suggestion unchanged, override it, or reject it entirely, is recorded in the subject log.

QC outputs are designed to be reviewed before making decisions. The pipeline pauses and waits for input at each manual prompt.

---

## QC Output Locations
| Location | Contents |
|----------|----------|
`sub-XXX/LOGS/` | All QC plots, metrics CSVs, channel map, labels map, main log |
`sub-XXX/QC/sub-XXX_icqc` | One IC QC packet figure per suggested IC

## LOGS: Per-Subject Log File
**File**: `LOGS/sub-XXX_preproc.log`

Timestamped plain-text log covering the full preprocessing run. Every log entry has the format:
```text
[YYYY-MM-DD HH:MM:SS] [STAGE] message
```
The log records:
- Experiment ID, output prefix, output root
- Raw file path resolved for this subject
- RNG seed
- All stage parameters (filter cutoffs, notch frequency, reref mode, etc.)
- Whether each stage was loaded from cache or processed fresh
- Automated bad channel suggestions with reasons
- Manual bad channel interpolation decision (including channels suggested but not interpolated)
- ICA training copy details (segments removed, % time, interval count)
- ICA method, extended mode flag, trained on FULL or CLEAN-COPY
- ICLabel suggestions with probabilities
- Manual IC rejection decision (including ICs suggested but kept)
- Every type counts before epoching
- Warnings for any non-fatal failures (ELP missing, ICLabel plugin failure, etc.)
- `[SKIP]` entries when a cached stage is loaded

The log is the primary audit trail for every preprocessing decision.

---

## LOGS: Channel Map and Labels
`sub-XXX_channelmap_applied.tsv`

Written by the montage step when `cfg.exp.montage.enabled = true`.  
Records the final channel mapping applied.
```text
index   label
1       Fp1
2       Fp2
```

`sub-XXX_chanlabels.txt`

Written by INITREJ step. Maps channel index to label.
```text
1   Fp1
2   Fp2
```

Useful for interpreting the bar charts and interpolation prompts by index.

---

## LOGS: Channel PSD Metrics
**File**: `LOGS/sub-XXX_chan_psd_metrics.csv`</br>
Written during INITREJ. One row per channel.

| Column | Description |
|--------|-------------|
`chan_idx` | Channel index
`label` | Channel label
`line_ratio` | Power 59–61 Hz / power (55–59 + 61–65 Hz). High values indicate line noise
`hf_ratio` | Power 20–40 Hz / power 1–12 Hz. High values indicate muscle artifact
`drift_ratio` | Power 1–2 Hz / power 1–12 Hz. High values indicate slow drift
`alpha_ratio` | Power 8–12 Hz / power 1–40 Hz. Proportion of broadband power in alpha band

All ratios are computed via `trapz` (trapezoidal integration of the Welch PSD). Welch parameters: 2-second windows, 50% overlap, FFT size = next power of 2 above window length.

---

## LOGS: INITREJ QC Plots
Seven figure types are saved during INITREJ:

### Histograms
**Files**: `sub-XXX_hist_chan_std.png`, `sub-XXX_hist_chan_rms.png`</br>
Distribution of per-channel STD and RMS across all channels. Outliers should be obvious.

### Bar Charts
**Files**: `sub-XXX_initrej_bar_chan_std.png`, `sub-XXX_initrej_bar_chan_rms.png`</br>
Per-channel STD and RMS as bar charts with suggested bad channels marked by vertical dashed lines.

### Topoplots
**Files**: `sub-XXX_initrej_topo_std.png`, `sub-XXX_initrej_topo_rms.png`</br>
Spatial distribution of STD and RMS on the scalp. Requires channel coordinates (ELP lookup). Useful for identifying topographically clustered artifacts.

### Channel PSD Topoplots
**Files**: `sub-XXX_initej_topo_line_ratio.png`, `sub-XXX_initej_topo_hf_ratio.png`, `sub-XXX_initej_topo_drift_ratio.png`, `sub-XXX_initej_topo_alpha_ratio.png`</br>
Topoplots of each PSD metric across channels. Line noise and HF ratio maps are particularly useful for identifying artifact-contaminated channels with normal amplitude statistics.

### PSD Overview
**File**: `sub-XXX_initrej_psd_overview.png`</br>
Median <u>+</u> IQR (25th/75th percentile) PSD across all channels (0–80 Hz). Provides a quick sanity check on the spectral profile of the cleaned data.

### Bad Channel PSD Overlay
**File**: `sub-XXX_initrej_psd_badchans.png`</br>
PSD for each suggested bad channel overlaid on a single plot. Only present when bad channels are suggested. Useful for confirming whether suggestions are genuine artifacts (flat line, extreme broadband elevation, line noise spike).

---

## INITREJ: Manual Bad Channel Decision
The pipeline prints suggested bad channels with their reasons and waits for input:
```text
[INITREJ]: Suggested bad channels: [3 17 42]
    1) 3 (Fp1): probability > 5 (pop_rejchan)
    2) 17 (F7): kurtosis > 5 (pop_rejchan)
    3) 42 (T7): probability > 5 + kurtosis > 5 (pop_rejchan)

Tyoe channel indices to INTERPOLATE (e.g., [1 2 17]).
Default = none (Press Enter or type []).
Channels to interpolate:
```
**Decision Options**:
- Enter the suggested indices to accept recommendation
- Enter a different set of indices to override
- Press Enter or type `[]` to interpolate nothing
- Enter additional indices beyond those suggested

The pipeline logs all three cases:
- Channels interpolated
- Channels suggested but not interpolated
- Interpolation method: spherical (`pop_interp`)

Channels are interpolated spherically from surrounding channels. The EEG retains its original channel count, interpolated channels remain in the data but are reconstructed.

---

## QC: IC QC Packets
**Location**: `QC/sub-XXX_icqc/sub-XXX_icNNN_qc.png`</br>
One figure per IC suggested for rejection by ICLabel. Each figure has three panels:

**Top left – Scalp Map**: Topoplot of the IC spatial filter (`EEG.icawinv(:, ic)`). Eye ICs show frontal polarity; muscle ICs show peripheral high-amplitude foci.

**Top right – Power Spectrum**: Welch PSD of the IC activation (0–80 Hz, dB scale). Muscle ICs show broad HF elevation; eye ICs show elevated low-frequency power; line noise ICs show a 60 Hz spike.

**Bottom – Activation time series**: First 10 seconds of the IC activation signal. Blink ICs show sterotyped slow deflections; muscle shows high-frequency bursts.

**Title**: ICLabel probability scores for all 7 categories; Brain (B), Muscle (M), Eye (E), Heart (H), Line Noise (L), Channel Noise (C), Other (O).

These figures should be opened and reviewed **before** entering the IC rejection prompt.

---

## LOGS: IC PSD Metrics
**File**: `LOGS/sub-XXX_ic_psd_metrics.csv`</br>
Written for suggested ICs after ICLabel. One row per suggested IC.

Column | Description
-------|------------
`ic` | IC index
`peak_hz` | Frequency of peak power (0.5–40 Hz)
`delta` | Delta band power (1–4 Hz) via trapz
`theta` | Theta band power (4–8 Hz)
`alpha` | Alpha band power (8–12 Hz)
`beta` | Beta band power (13–30 Hz)
`gamma` | Gamma band power (30–45 Hz)
`hf_ratio` | Power 20–40 Hz / power 1–12 Hz
`line_ratio` | Line noise ratio (59–61 vs. flanking bands)

---

## ICA: Manual IC Rejection Decision
After reviewing QC packets, the pipeline prompts:
```text
[ICREJ] ICLabel suggested ICs: [1 3 7]
Review QC figs in QC/sub-XXX_icqc/ before deciding.
Type IC indices to REMOVE (e.g., [1 3 7]).
Default = remove none (press Enter or type []).
ICs to remove:
```
- The default is to remove **nothing**; conservative by design
- Entering indices removes them via `pop_subcomp`
- ICs suggested but kept are listed in the log with their ICLabel reasons

---

## ICA Training Copy: Bad Segment Prompt
If `badseg.enabled = true`, the pipeline detects segments exceeding a voltage threshold and prompts:
```text
[ICA-TRAIN] Detected 12 intervals (3.24% of samples).
Remove detected bad segments from ICA training copy? (y/n) [n]:
```
- Default is **no**; keeps all data
- If accepted, segments are removed only from the ICA training copy; the full length EEG is preserved
- The number of intervals, percentage of time, and whether removal occured are all logged

---

## Interpreting the Log
Key patterns to search for in `sub-XXX_preproc.log`:

Pattern | Meaning
--------|---------
`[SKIP]` | Stage loaded from cache; not reprocessed
`[WARN]` | Non-fatal warning; pipeline continued
`[INITREJ] Suggested bad channels: []` | No channels flagged by automated methods
`[INITREJ] No channels interpolated` | Manual decision to skip interpolation
`[ICREJ] No ICs removed` | Manual decision to keep all suggested ICs
`[ICREJ] Suggested but NOT removed:` | ICLabel flagged but human disagreed
`[INITREJ] Suggested but NOT interpolated` | Pipeline flagged but human disagreed
`[ICA-TRAIN] Trained on CLEAN-COPY` | ICA was trained on data with bad segments removed