# Filesystem
This page documents the complete directory structure produced by the framework, the naming conventions for every file type, and the separation between read-only raw data and writable variables.

---

## Root Structure
```
/cifs/seminowicz/eegPainDatasets/CNED/                 ← RAW_ROOT (read-only)
│   └── <raw_dirname>/                                  e.g. 26ByBiosemi/
│       ├── sub-001/eeg/sub-001_task-*_eeg.bdf
│       ├── sub-002/eeg/...
│       └── participants.tsv
│
└── da-analysis/                                        ← PROJ_ROOT (writable)
    └── <out_dirname>/                                  e.g. 26ByBiosemi/
        ├── resource/                                   ← P.RESOURCE
        │   ├── standard-10-5-cap385.elp
        │   ├── montage.csv
        │   └── participants.tsv
        └── preproc/                                    ← P.RUN_ROOT
            └── sub-001/                                per-subject root
                ├── 01_filter/
                ├── 02_notch/
                ├── 03_resample/
                ├── 04_reref/
                ├── 05_initrej/
                ├── 06_ica/
                ├── 07_epoch/
                ├── 08_base/
                ├── LOGS/
                ├── QC/
                └── SPECTRAL/                           (created by spectral pipeline)
                    ├── csv/
                    ├── figures/
                    ├── tmp/
                    └── logs/
```

Raw data is never written to. All outputs live under <code>PROJ_ROOT</code>.

---

## Per-Subject Stage Directories
Each stage directly holds <code>.set</code> files for that subject at that preprocessing step. The filename encodes the full preprocessing history via culmulative tags.
```
01_filter/
    26BB_62_001_fir.set
    26BB_62_001_fir.fdt

02_notch/
    26BB_62_001_fir_notch60.set
    26BB_62_001_fir_notch60.fdt

03_resample/
    26BB_62_001_fir_notch60_rs500.set       (only if resample enabled)

04_reref/
    26BB_62_001_fir_notch60_reref.set

05_initrej/
    26BB_62_001_fir_notch60_reref_initrej.set

06_ica/
    26BB_62_001_fir_notch60_reref_initrej_ica.set
    26BB_62_001_fir_notch60_reref_initrej_ica_iclabel.set   (if iclabel tag set)

07_epoch/
    26BB_62_001_fir_notch60_reref_initrej_ica_epoch.set

08_base/
    26BB_62_001_fir_notch60_reref_initrej_ica_epoch_base.set
```

---

## Filename Convention
```
<out_prefix><subjid_zero_padded_3>_<tag1>_<tag2>_..._<tagN>.set
```

The naming function is <code>P.NAMING.fname(subjid, tags, prefix)</code> -> <code>local_fname.m</code>.

---

## LOGS Directory

One LOGS directory per subject: <code>sub-XXX/LOGS/</code>.
```
LOGS/
    sub-001_preproc.log                 Main timestamped log (all stages)

    sub-001_chalabels.txt               Channel index → label map
    sub-001_channelmap_applied.tsv      Montage audit (raw_label → std_label)

    sub-001_chan_psd_metrics.csv        Per-channel PSD metrics (INITREJ)
    sub-001_ic_psd_metrics.csv          Per-IC PSD metrics (ICA/ICLabel)

    sub-001_initrej_hist_chan_std.png    STD histogram (all channels)
    sub-001_initrej_hist_chan_rms.png    RMS histogram (all channels)
    sub-001_initrej_bar_chan_std.png     STD bar chart (bad chans marked)
    sub-001_initrej_bar_chan_rms.png     RMS bar chart (bad chans marked)
    sub-001_initrej_topo_std.png        STD topoplot
    sub-001_initrej_topo_rms.png        RMS topoplot
    sub-001_initrej_psd_overview.png    Median ± IQR PSD (all channels)
    sub-001_initrej_psd_badchans.png    PSD overlay (suggested bad channels only)

    sub-001_initrej_topo_line_ratio.png   Line noise topoplots
    sub-001_initrej_topo_hf_ratio.png
    sub-001_initrej_topo_drift_ratio.png
    sub-001_initrej_topo_alpha_ratio.png
```

---

## QC Directory
One QC directory per subject: <code>sub-XXX/QC/</code>

```
QC/
    sub-001_icqc/
        sub-001_ic001_qc.png    IC QC packet (scalp map + PSD + activation)
        sub-001_ic002_qc.png
        ...                     One figure per suggested IC
```

QC packets are generated for all ICs suggested by ICLabel. They are reviewed before making manual IC rejection decisions.

---

## SPECTRAL Directory
Created by the spectral pipeline: <code>sub-XXX/SPECTRAL</code>.

```
SPECTRAL/
    csv/
        sub-001_spectral_chan_by_trial.csv    Per-channel × per-trial alpha features
        sub-001_spectral_ga_by_trial.csv      GA (averaged across channels) × per-trial

    figures/
        sub-001_spectral_summary.png          GA PSD + PAF + alpha interaction summary
        sub-001_spectral_fooof.png            FOOOF trial-wise summary
        sub-001_heatmap_panel.png             Channel × trial heatmaps (all features)

    trial_spectral/                           (if trial_spectral.enabled)
        sub-001_trial-001_prepost_psd.png     Pre vs post-stim PSD per trial
        sub-001_trial-002_prepost_psd.png
        ...

    tmp/
        sub-001_fooof_freqs.csv     FOOOF input: frequency vector
        sub-001_fooof_psd.csv       FOOOF input: GA PSD per trial
        sub-001_fooof_cfg.json      FOOOF parameters passed to Python
        sub-001_fooof_out.json      Raw FOOOF output from Python

    logs/
        sub-001_spectral.log        Spectral pipeline log (timestamped)
```

---

## Resource Directory
<code><out_dirname>/resource/</code> holds files shared across all subjects in an experiment.

```
resource/
    standard-10-5-cap385.elp       Channel location coordinates (ELP format)
    montage.csv                    Biosemi A/B label → standard label mapping
    participants.tsv               Subject list and demographics
    cap_size.csv                   EEG cap size per experiment
```

The ELP file is exposed as <code>P.CORE.ELP_FILE</code>. The participants TSV is <code>P.CORE.PARTICIPANTS_TSV</code>.

---

## R Analysis Outputs
```
R/analysis/behavioural_analysis/
    behavioural_master.csv                  All experiments, all trials, raw behaviour
    behavioural_demo_master.csv             + age, sex, cap_size merged in
    alpha_pain_master.csv                   + spectral EEG features merged in

    alpha_pain_master_model_input.csv       Filtered + z-scored model input
    alpha_pain_master_model_input_v2.csv    v2 model input

    model_comparison.csv                    AIC/BIC/logLik for all models
    model_comparison_v2.csv

    trial_level_fitted_values.csv           Per-trial fitted + residuals (all models)
    trial_level_fitted_values_v2.csv
    subject_level_ga_fitted_values.csv      Subject-level grand-average fitted values
    subject_level_ga_fitted_values_v2.csv

    gamm_outputs/
        mXX_baseline_summary.txt            summary(mod) for each model
        mXX_baseline_diagnostics.png        gam.check() plots
        mXX_baseline_observed_vs_fitted.png
        mXX_baseline.rds                    Serialized model object

    gamm_outputs_v2/
        ...                                 Mirror structure for v2 models
        m05_tensor_slow_fast_surface.png    Slow×fast interaction contour
        m09_tensor_slow_fast_aperiodic_surface.png
```

---

## Config Directory
```
config/
    exp01.json
    exp02.json
    ...
```

---

## Path Override via JSON
All root paths can be overidden in the JSON config without editing any MATLAB code:
```json
{
    "paths": {
        "raw_root": "/my/custom/raw/path",
        "proj_root": "/my/custom/output/path"
    }
}
```

<code>config_paths.m</code> reads these via <code>cfg_get_string()</code> with sensible defaults.