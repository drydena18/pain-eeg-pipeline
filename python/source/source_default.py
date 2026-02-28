from __future__ import annotations

import os
import json
from copy import deepcopy

from source_core import source_core

def _get(d, path, default = None):
    cur = d
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur
    
def _ensure_dir(p: str):
    os.makedirs(p, exist_ok = True)

def source_default(exp_id: str, cfg_in: dict, subjects_override = None):
    cfg = deepcopy(cfg_in)

    if "exp" not in cfg:
        raise ValueError("Missing cfg.exp in JSON")
    if "source" not in cfg:
        raise ValueError("Missing cfg.source in JSON")
    
    # exp_id fallback
    cfg["exp"]["id"] = str(cfg["exp"].get("id", exp_id))

    # subjects
    if subjects_override is not None and len(subjects_override) > 0:
        cfg["exp"]["subjects"] = [int(s) for s in subjects_override]
    else:
        if "subjects" not in cfg["exp"] or not cfg["exp"]["subjects"]:
            raise ValueError("cfg.exp.subjects is empty. (For now) please provide subjects explicitly.")
        cfg["exp"]["subjects"] = [int(s) for s in cfg["exp"]["subjects"]]

    # ---- Source defaults ----
    src = cfg["source"]
    src["enabled"] = bool(src.get("enabled", True))

    # input
    src.setdefault("input", {})
    src["input"].setdefault("stage_dir", "08_base")
    src["input"].setdefault("allow_fallback_search", True)

    # paths: da-anaysis root
    # Let JSON define da_root and exp_out
    da_root = _get(cfg, ["paths", "da_root"], "/cifs/seminowic/eegPainDatasets/CNED/da-analysis")
    exp_out = _get(cfg, ["exp", "out_dirname"], None) or _get(cfg, ["exp", "id"], None)
    if exp_out is None:
        raise ValueError("Need cfg.exp.out_dirname (recommended) or cfg.exp.id to infer exp_out folder name.")
    
    # output root
    src.setdefault("outputs", {})
    out_root = src["outputs"].get("root", "AUTO")
    if out_root == "AUTO":
        out_root = os.path.join(da_root, exp_out, "source")
    src["outputs"]["root"] = out_root

    # fsaverage assets (bundled inside repo or on server)
    src.setdefault("fsaverage", {})
    if "subjects_dir" not in src["fsaverage"] or not src["fsaverage"]["subjects_dir"]:
        raise ValueError('cfg.source.fsaverage.subjects_dir is required (path that contains "fsaverage/").')
    
    # Core MNE params
    src.setdefault("lambda2", 1.0 / 9.0)
    src.setdefault("forward", {})
    src["forward"].setdefault("mindist_mm", 5.0)

    src.setdefault("inverse", {})
    src["inverse"].setdefault("method", "sLORETA")
    src["inverse"].setdefault("snr", 3.0)
    src["inverse"].setdefault("pick_ori", None)

    src.setdefault("noise_cov", {})
    src["noise_cov"].setdefault("tmin", -0.2)
    src["noise_cov"].setdefault("tmax", 0.0)

    # ROI defaults
    src.setdefault("roi", {})
    src["roi"].setdefault("parcellation", "aparc")
    src["roi"].setdefault("mode", "mean_flip")
    src["roi"].setdefault("use_custom_rois", True)
    src["roi"].setdefault("custom_rois", [])

    # Spectral defaults for source-space features
    src.setdefault("spectral", {})
    src["spectral"].setdefault("fmin", 1.0)
    src["spectral"].setdefault("fmax", 40.0)
    src["spectral"].setdefault("alpha_band", [8.0, 12.0])
    src["spectral"].setdefault("slow_alpha_band", [8.0, 10.0])
    src["spectral"].setdefault("fast_alpha_band", [10.0, 12.0])

    # FOOOF defaults
    src.setdefault("fooof", {})
    src["fooof"].setdefault("enabled", True)
    src["fooof"].setdefault("mode", "ga-only") # "ga-only" or "trial-and-ga"
    src["fooof"].setdefault("aperiodic_mode", "fixed") # fixed or knee
    src["fooof"].setdefault("peak_width_limits", [1.0, 12.0])
    src["fooof"].setdefault("max_n_peaks", 6)
    src["fooof"].setdefault("min_peak_height", 0.1)
    src["fooof"].setdefault("peak_threshold", 2.0)
    src["fooof"].setdefault("freq_range", [1.0, 40.0])

    # QC defaults
    src.setdefault("qc", {})
    src["qc"].setdefault("save_plots", True)
    src["qc"].setdefault("save_brain_images", True)
    src["qc"].setdefault("brain_snapshot_time_sec", 0.05)

    _ensure_dir(src["outputs"]["root"])

    print(f"[{cfg['exp']['id']}] SOURCE subjects ({len(cfg['exp']['subjects'])}): {cfg['exp']['subjects']}")
    print(f"Source input stage_dir: {src['input']['stage_dir']}")
    print(f"Source output root: {src['outputs']['root']}")

    if src["enabled"]:
        source_core(cfg, da_root = da_root, exp_out = exp_out)
    else:
        print(f"[{cfg['exp']['id']}] cfg.source.enabled = false (skipping)")