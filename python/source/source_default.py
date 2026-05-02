"""
source_default.py  –  Config validation, safe defaults, and dispatch.
V 2.0.0

Call chain:
    exp01_source.py  ->  source_default()  ->  source_core()

Mirrors spectral_default.m:
    - Hard errors for genuinely required fields
    - Safe defaults for everything else (no JSON changes needed for a basic run)
    - Run-header print before dispatch

Config key: cfg["source"]  (parallel to cfg["spectral"])
"""

from __future__ import annotations

import os
from copy import deepcopy

from source_core import source_core


# =============================================================================
# HELPERS
# =============================================================================

def _d(d: dict, key: str, val):
    """Set d[key] = val only if key is absent (equivalent to defaultField)."""
    d.setdefault(key, val)
    return d


def _require(d: dict, key: str, label: str):
    """Raise a clear ValueError if a required key is missing or empty."""
    v = d.get(key, None)
    if v is None or v == "" or v == [] or v == {}:
        raise ValueError(f"Required config field missing or empty: {label}")


def _get(d: dict, path: list, default=None):
    """Safe nested dict access."""
    cur = d
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


# =============================================================================
# ENTRY POINT
# =============================================================================

def source_default(exp_id: str, cfg_in: dict, subjects_override=None):
    """
    Validate config, fill safe defaults, then dispatch to source_core().

    Args:
        exp_id            : Experiment identifier (e.g. "exp01").
        cfg_in            : Raw dict loaded from the experiment JSON.
        subjects_override : Optional iterable of integer subject IDs.
                            Overrides cfg["exp"]["subjects"] when provided.
    """
    cfg = deepcopy(cfg_in)

    # ── Required top-level blocks ─────────────────────────────────────────────
    _require(cfg, "exp",    "cfg.exp")
    _require(cfg, "source", "cfg.source")

    # ── Experiment block ──────────────────────────────────────────────────────
    exp = cfg["exp"]
    exp["id"] = str(exp.get("id", exp_id))

    _require(exp, "out_prefix",
             "cfg.exp.out_prefix  (e.g. '26BB_62_' — needed to resolve .set filenames)")

    if subjects_override is not None and len(list(subjects_override)) > 0:
        exp["subjects"] = [int(s) for s in subjects_override]
    else:
        _require(exp, "subjects", "cfg.exp.subjects (list of integer subject IDs)")
        exp["subjects"] = [int(s) for s in exp["subjects"]]

    if len(exp["subjects"]) == 0:
        raise ValueError("cfg.exp.subjects resolved to an empty list.")

    # ── Path resolution ───────────────────────────────────────────────────────
    da_root = _get(cfg, ["paths", "da_root"],
                   "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis")
    exp_out = exp.get("out_dirname") or exp.get("id")
    if not exp_out:
        raise ValueError(
            "Cannot resolve experiment output folder. "
            "Provide cfg.exp.out_dirname (recommended) or cfg.exp.id."
        )

    # ── Source block defaults ─────────────────────────────────────────────────
    src = cfg["source"]
    _d(src, "enabled", True)

    # Input
    _d(src, "input", {})
    _d(src["input"], "stage_dir",            "08_base")
    _d(src["input"], "allow_fallback_search", True)

    # Output root
    _d(src, "outputs", {})
    if src["outputs"].get("root", "AUTO") == "AUTO":
        src["outputs"]["root"] = os.path.join(da_root, exp_out, "source")

    # fsaverage assets (no safe default — must be explicit)
    _d(src, "fsaverage", {})
    _require(src["fsaverage"], "subjects_dir",
             "cfg.source.fsaverage.subjects_dir  "
             "(path containing fsaverage/; obtain via mne.datasets.fetch_fsaverage())")

    # Forward solution
    _d(src, "forward", {})
    _d(src["forward"], "mindist_mm", 5.0)

    # Inverse operator
    _d(src, "inverse", {})
    _d(src["inverse"], "method",   "sLORETA")
    _d(src["inverse"], "snr",      3.0)
    _d(src["inverse"], "loose",    0.2)
    _d(src["inverse"], "depth",    0.8)
    _d(src["inverse"], "pick_ori", None)

    # Noise covariance baseline window
    _d(src, "noise_cov", {})
    _d(src["noise_cov"], "tmin", -0.2)
    _d(src["noise_cov"], "tmax",  0.0)

    # Pre-stimulus window
    _d(src, "prestim", {})
    _d(src["prestim"], "tmin", -0.5)
    _d(src["prestim"], "tmax",  0.0)

    # Post-stimulus window
    _d(src, "poststim", {})
    _d(src["poststim"], "tmin",        0.1)   # avoid stimulus artefact
    _d(src["poststim"], "tmax",        0.8)
    _d(src["poststim"], "phase_ref_t", 0.2)   # time at which post-stim phase is sampled

    # LEP (laser-evoked potential) windows
    _d(src, "lep", {})
    _d(src["lep"], "n2_window", [0.15, 0.35])
    _d(src["lep"], "p2_window", [0.25, 0.50])

    # ROI (region of interest) settings
    # Default custom ROIs cover the core pain neuromatrix in the Desikan-Killiany
    # (aparc) parcellation.  Labels are bilateral (lh + rh merged).
    # To use all parcels instead, set use_custom_rois = false.
    #
    # Parcellation label reference (aparc):
    #   S1  — postcentral         (primary somatosensory cortex)
    #   M1  — precentral          (primary motor cortex)
    #   ACC — caudalanteriorcingulate + rostralanteriorcingulate
    #   Ins — insula               (insular cortex)
    #   SII — supramarginal        (secondary somatosensory / parietal operculum proxy)
    #   dlPFC — rostralmiddlefrontal  (dorsolateral prefrontal — descending modulation)
    _d(src, "roi", {})
    _d(src["roi"], "parcellation",    "aparc")
    _d(src["roi"], "mode",            "mean_flip")
    _d(src["roi"], "use_custom_rois",  True)
    _d(src["roi"], "custom_rois", {
        "S1":   ["postcentral-lh",              "postcentral-rh"],
        "M1":   ["precentral-lh",               "precentral-rh"],
        "ACC":  ["caudalanteriorcingulate-lh",   "caudalanteriorcingulate-rh",
                 "rostralanteriorcingulate-lh",  "rostralanteriorcingulate-rh"],
        "Ins":  ["insula-lh",                   "insula-rh"],
        "SII":  ["supramarginal-lh",             "supramarginal-rh"],
        "dlPFC":["rostralmiddlefrontal-lh",      "rostralmiddlefrontal-rh"],
    })

    # Spectral parameters
    _d(src, "spectral", {})
    _d(src["spectral"], "fmin",             1.0)
    _d(src["spectral"], "fmax",            40.0)
    _d(src["spectral"], "alpha_band",      [8.0, 12.0])
    _d(src["spectral"], "slow_alpha_band", [8.0, 10.0])
    _d(src["spectral"], "fast_alpha_band", [10.0, 12.0])

    # FOOOF (Fitting Oscillations & One Over F)
    _d(src, "fooof", {})
    _d(src["fooof"], "enabled",           True)
    _d(src["fooof"], "aperiodic_mode",    "fixed")   # "fixed" or "knee"
    _d(src["fooof"], "peak_width_limits", [1.0, 12.0])
    _d(src["fooof"], "max_n_peaks",       6)
    _d(src["fooof"], "min_peak_height",   0.1)
    _d(src["fooof"], "peak_threshold",    2.0)
    _d(src["fooof"], "freq_range",        [1.0, 40.0])

    # Quality control / 2-D plots
    _d(src, "qc", {})
    _d(src["qc"], "save_plots", True)

    # 3-D brain rendering
    _d(src, "render", {})
    _d(src["render"], "enabled",      False)   # opt-in; requires pyvista
    _d(src["render"], "use_mesa",     False)   # True on headless HPC servers
    _d(src["render"], "stc_enabled",  True)    # render GA STC at pre/N2/P2
    _d(src["render"], "stc_clim_pct", [50, 99])  # percentile colour scale bounds
    _d(src["render"], "roi_enabled",  True)    # render ROI scalar maps
    _d(src["render"], "roi_metrics",          # which GA metrics to paint
       ["BI_pre", "LR_pre", "CoG_pre", "delta_ERD", "n2p2_amp"])

    # Validate inverse method
    if str(src["inverse"]["method"]).lower() != "sloreta":
        raise ValueError(
            f"cfg.source.inverse.method = '{src['inverse']['method']}' is not supported. "
            "Only 'sLORETA' is currently implemented."
        )

    os.makedirs(src["outputs"]["root"], exist_ok=True)

    # ── Run header ────────────────────────────────────────────────────────────
    print(f"\n[{exp['id']}] SOURCE  subjects ({len(exp['subjects'])}): {exp['subjects']}")
    print(f"  Input stage    : {src['input']['stage_dir']}")
    print(f"  Output root    : {src['outputs']['root']}")
    print(f"  Parcellation   : {src['roi']['parcellation']}"
          + (" (custom ROIs)" if src["roi"]["use_custom_rois"] else ""))
    print(f"  Pre-stim       : {src['prestim']['tmin']:.3f} – {src['prestim']['tmax']:.3f} s")
    print(f"  Post-stim      : {src['poststim']['tmin']:.3f} – {src['poststim']['tmax']:.3f} s "
          f"  (phase ref @ {src['poststim']['phase_ref_t']:.3f} s)")
    print(f"  Noise cov      : {src['noise_cov']['tmin']:.3f} – {src['noise_cov']['tmax']:.3f} s")
    print(f"  LEP N2 window  : {src['lep']['n2_window']}")
    print(f"  LEP P2 window  : {src['lep']['p2_window']}")
    print(f"  FOOOF          : {'enabled (' + src['fooof']['aperiodic_mode'] + ')' if src['fooof']['enabled'] else 'disabled'}")
    if src["render"]["enabled"]:
        print(f"  3-D renders    : STC={'on' if src['render']['stc_enabled'] else 'off'}  "
              f"ROI={'on' if src['render']['roi_enabled'] else 'off'}  "
              f"Mesa={'on' if src['render']['use_mesa'] else 'off'}  "
              f"metrics={src['render']['roi_metrics']}")

    # ── Dispatch ──────────────────────────────────────────────────────────────
    if src["enabled"]:
        source_core(cfg, da_root=da_root, exp_out=exp_out)
    else:
        print(f"[{exp['id']}] cfg.source.enabled = false — skipping.")