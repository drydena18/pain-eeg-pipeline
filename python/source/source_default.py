"""
source_default.py  -  Config validation, safe defaults, and dispatch.
V 3.1.0

Call chain:
    expXX_source.py  ->  source_default()  ->  source_core()

V3.1.0: features.windows.whole_sec (default [-1, 2] s) replaces "whole =
full epoch" now that epochs run -2.2 to 2.0 s; validation warns when a window
lies within half a band-power filter length of an epoch edge.

V3.0.0 changes vs V2.1.0:
    - Feature parameters now come from the top-level cfg["features"] block,
      the single source of truth shared with spectral_default.m (windows,
      bands, band-power filter, Welch PSD, ERD guard, LEP windows, phase
      reference, FOOOF). They are mapped onto the source.* keys that
      source_core reads. A conflicting legacy source.* value triggers a
      warning and the features value wins.
    - Features are validated against the PREPROCESSING baseline in the same
      JSON (epoch limits, band-pass filter): windows must lie inside the
      epoch, bands inside the pass band, and the band-power FIR must fit in
      the epoch. Identical checks run in spectral_default.m.
    - ROI defaults aligned with the experiment configs (S1 includes
      paracentral, S2 = supramarginal, Insula, dlPFC = rostral + caudal
      middle frontal); roi.hemisphere_split (default true) adds _lh / _rh
      versions of every bilateral ROI.
    - Removed quiet_band (see src_erd.py V3.0.0); added eps0_frac / p5_percentile.
"""

from __future__ import annotations

import math
import os
from copy import deepcopy

from source_core import source_core


# =============================================================================
# SHARED FEATURE DEFAULTS (mirror of spectral_default.m -> local_feature_defaults)
# =============================================================================
FEATURE_DEFAULTS = {
    "windows":    {"whole_sec": [-1.0, 2.0], "pre_sec": [-1.0, -0.1], "post_sec": [0.1, 0.8]},
    "bands":      {"alpha_hz": [8.0, 12.0], "slow_hz": [8.0, 10.0], "fast_hz": [10.0, 12.0]},
    "band_power": {"method": "filter_hilbert", "trans_bw_hz": 1.5},
    "psd":        {"fmin_hz": 1.0, "fmax_hz": 40.0, "window_sec": 2.0,
                   "overlap_frac": 0.5, "df_target_hz": 0.25},
    "erd":        {"eps0_frac": 1e-3, "p5_percentile": 5.0},
    "lep":        {"n2_window_sec": [0.15, 0.35], "p2_window_sec": [0.25, 0.50]},
    "phase":      {"post_ref_sec": 0.2},
    "fooof":      {"enabled": True, "aperiodic_mode": "fixed", "peak_width_limits": [1.0, 12.0],
                   "max_n_peaks": 6, "min_peak_height": 0.1, "peak_threshold": 2.0,
                   "freq_range": [1.0, 40.0]},
}

DEFAULT_CUSTOM_ROIS = {
    "S1":     ["lh-postcentral", "rh-postcentral", "lh-paracentral", "rh-paracentral"],
    "S2":     ["lh-supramarginal", "rh-supramarginal"],
    "ACC":    ["lh-caudalanteriorcingulate", "rh-caudalanteriorcingulate",
               "lh-rostralanteriorcingulate", "rh-rostralanteriorcingulate"],
    "Insula": ["lh-insula", "rh-insula"],
    "dlPFC":  ["lh-rostralmiddlefrontal", "rh-rostralmiddlefrontal",
               "lh-caudalmiddlefrontal", "rh-caudalmiddlefrontal"],
    "M1":     ["lh-precentral", "rh-precentral"],
}


# =============================================================================
# HELPERS
# =============================================================================
def _d(d: dict, key: str, val):
    """Set d[key] = val only if key is absent (equivalent to defaultField)."""
    d.setdefault(key, val)
    return d


def _require(d: dict, key: str, label: str):
    v = d.get(key, None)
    if v is None or v == "" or v == [] or v == {}:
        raise ValueError(f"Required config field missing or empty: {label}")


def _get(d: dict, path: list, default = None):
    cur = d
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def _merge_defaults(user: dict, defaults: dict) -> dict:
    out = deepcopy(defaults)
    for k, v in (user or {}).items():
        if k.startswith("_"):
            continue
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge_defaults(v, out[k])
        else:
            out[k] = v
    return out


def _set_from_features(block: dict, key: str, value, label: str):
    """Write a features value into a legacy source.* key, warning on conflict."""
    if key in block and block[key] != value:
        print(f"[WARN] {label} = {block[key]} conflicts with cfg.features ({value}); using features.")
    block[key] = value


def resolve_features(cfg: dict) -> dict:
    """Return cfg['features'] filled with defaults (no validation)."""
    return _merge_defaults(cfg.get("features", {}), FEATURE_DEFAULTS)


def validate_features_against_preproc(feat: dict, cfg: dict, sfreq: float | None = None):
    """
    Check the shared features against the preprocessing baseline.
    Raises ValueError on hard conflicts, prints warnings on soft ones.
    Mirrors spectral_default.m -> local_validate_features.
    """
    pp = cfg.get("preproc", {})
    ep = pp.get("epoch", {})
    tmin, tmax = float(ep.get("tmin_sec", -1.0)), float(ep.get("tmax_sec", 2.0))
    filt = pp.get("filter", {})
    hp, lp = float(filt.get("highpass_hz", 0.5)), float(filt.get("lowpass_hz", 40.0))
    fs = float(sfreq or _get(pp, ["resample", "target_hz"], 500) or 500)

    w = feat["windows"]
    for name in ("whole_sec", "pre_sec", "post_sec"):
        a, b = w[name]
        if not (tmin <= a < b <= tmax):
            raise ValueError(f"features.windows.{name} = {w[name]} is not inside the epoch [{tmin}, {tmax}] s.")

    b = feat["bands"]
    if not (b["slow_hz"][1] == b["fast_hz"][0] and b["slow_hz"][0] == b["alpha_hz"][0]
            and b["fast_hz"][1] == b["alpha_hz"][1]):
        print(f"[WARN] features.bands: slow {b['slow_hz']} + fast {b['fast_hz']} do not tile alpha {b['alpha_hz']}.")
    for name in ("alpha_hz", "slow_hz", "fast_hz"):
        lo, hi = b[name]
        if not (hp < lo < hi < lp):
            raise ValueError(f"features.bands.{name} = {b[name]} is outside the preprocessing pass band ({hp}-{lp} Hz).")

    tb = float(feat["band_power"]["trans_bw_hz"])
    n_taps = int(math.floor(3.3 * fs / tb + 0.5)) | 1
    n_epoch = int(round((tmax - tmin) * fs)) + 1
    if n_taps > n_epoch:
        raise ValueError(f"Band-power FIR ({n_taps} taps at {fs} Hz) is longer than the epoch ({n_epoch} samples); "
                         "raise features.band_power.trans_bw_hz.")
    # Edge zone: samples within half a filter length of an epoch edge depend on
    # the mirror padding (pre-window power ~9 % low when the window touches the edge).
    half = (n_taps - 1) / 2.0 / fs
    for name in ("whole_sec", "pre_sec", "post_sec"):
        a, b = w[name]
        if a - tmin < half - 1e-9 or tmax - b < half - 1e-9:
            print(f"[WARN] features.windows.{name} = {w[name]} lies within {half:.2f} s of the epoch edge "
                  f"[{tmin}, {tmax}] s; band power there depends on filter padding.")

    if float(feat["psd"]["fmax_hz"]) > lp:
        print(f"[WARN] features.psd.fmax_hz = {feat['psd']['fmax_hz']} exceeds the {lp} Hz low-pass.")

    post = w["post_sec"]
    for name in ("n2_window_sec", "p2_window_sec"):
        lo, hi = feat["lep"][name]
        if not (post[0] <= lo < hi <= post[1]):
            print(f"[WARN] features.lep.{name} = {[lo, hi]} is not inside the post window {post}.")

    ref = float(feat["phase"]["post_ref_sec"])
    if not (tmin < ref < tmax) or math.isnan(ref):
        raise ValueError(f"features.phase.post_ref_sec = {ref} is outside the epoch.")


# =============================================================================
# ENTRY POINT
# =============================================================================
def source_default(exp_id: str, cfg_in: dict, subjects_override = None):
    """Validate config, fill safe defaults, then dispatch to source_core()."""
    cfg = deepcopy(cfg_in)

    _require(cfg, "exp",    "cfg.exp")
    _require(cfg, "source", "cfg.source")

    exp = cfg["exp"]
    exp["id"] = str(exp.get("id", exp_id))
    _require(exp, "out_prefix", "cfg.exp.out_prefix  (e.g. '26BB_62_')")

    if subjects_override is not None and len(list(subjects_override)) > 0:
        exp["subjects"] = [int(s) for s in subjects_override]
    else:
        _require(exp, "subjects", "cfg.exp.subjects")
        exp["subjects"] = [int(s) for s in exp["subjects"]]
    if len(exp["subjects"]) == 0:
        raise ValueError("cfg.exp.subjects resolved to an empty list.")

    da_root = _get(cfg, ["paths", "da_root"], "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis")
    exp_out = exp.get("out_dirname") or exp.get("id")
    if not exp_out:
        raise ValueError("Provide cfg.exp.out_dirname (recommended) or cfg.exp.id.")

    # ── Shared features (single source of truth) ─────────────────────────────
    feat = resolve_features(cfg)
    validate_features_against_preproc(feat, cfg)
    cfg["features"] = feat

    # ── Source block ─────────────────────────────────────────────────────────
    src = cfg["source"]
    _d(src, "enabled", True)

    _d(src, "input", {})
    _d(src["input"], "stage_dir",             "08_base")
    _d(src["input"], "allow_fallback_search", True)

    _d(src, "outputs", {})
    if src["outputs"].get("root", "AUTO") == "AUTO":
        src["outputs"]["root"] = os.path.join(da_root, exp_out, "source")

    _d(src, "fsaverage", {})
    fsa_dir = src["fsaverage"].get("subjects_dir", "")
    if not fsa_dir or fsa_dir == "NEEDS_PATH":
        raise ValueError("cfg.source.fsaverage.subjects_dir is not set (directory containing fsaverage/).")

    _d(src, "forward", {})
    _d(src["forward"], "mindist_mm", 5.0)

    _d(src, "inverse", {})
    _d(src["inverse"], "method",   "sLORETA")
    _d(src["inverse"], "snr",      3.0)
    _d(src["inverse"], "loose",    0.2)
    _d(src["inverse"], "depth",    0.8)
    _d(src["inverse"], "pick_ori", "normal")

    _d(src, "noise_cov", {})
    _d(src["noise_cov"], "tmin", -0.2)
    _d(src["noise_cov"], "tmax",  0.0)

    apply_features_to_source(src, feat)

    # ROIs
    _d(src, "roi", {})
    _d(src["roi"], "parcellation",     "aparc")
    _d(src["roi"], "mode",             "mean_flip")
    _d(src["roi"], "use_custom_rois",  True)
    _d(src["roi"], "hemisphere_split", True)
    _d(src["roi"], "custom_rois",      deepcopy(DEFAULT_CUSTOM_ROIS))

    _d(src, "qc", {})
    _d(src["qc"], "save_brain_images",       False)
    _d(src["qc"], "brain_snapshot_time_sec", 0.200)

    if str(src["inverse"]["method"]).lower() != "sloreta":
        raise ValueError(f"cfg.source.inverse.method = '{src['inverse']['method']}' is not supported (sLORETA only).")
    if src["forward"].get("spacing") not in (None, "ico5"):
        print(f"[WARN] source.forward.spacing = {src['forward']['spacing']} is ignored; "
              "the fsaverage ico-5 source space is always used.")

    os.makedirs(src["outputs"]["root"], exist_ok = True)

    # ── Run header ───────────────────────────────────────────────────────────
    roi = src["roi"]
    print(f"\n[{exp['id']}] SOURCE  subjects ({len(exp['subjects'])}): {exp['subjects']}")
    print(f"  Input stage    : {src['input']['stage_dir']}")
    print(f"  Output root    : {src['outputs']['root']}")
    print(f"  Parcellation   : {roi['parcellation']}"
          + (f" (custom ROIs{', + _lh/_rh' if roi['hemisphere_split'] else ''})" if roi["use_custom_rois"] else " (all labels)"))
    print(f"  Windows        : whole {feat['windows']['whole_sec']} s, pre {feat['windows']['pre_sec']} s, post {feat['windows']['post_sec']} s")
    print(f"  Band power     : filter-Hilbert, slow {feat['bands']['slow_hz']} / fast {feat['bands']['fast_hz']} Hz, "
          f"trans_bw {feat['band_power']['trans_bw_hz']} Hz")
    print(f"  Noise cov      : {src['noise_cov']['tmin']:.3f} - {src['noise_cov']['tmax']:.3f} s")
    print(f"  LEP N2 / P2    : {feat['lep']['n2_window_sec']} / {feat['lep']['p2_window_sec']} s")
    print(f"  FOOOF          : {'enabled (' + feat['fooof']['aperiodic_mode'] + ')' if feat['fooof']['enabled'] else 'disabled'}")

    if src["enabled"]:
        source_core(cfg, da_root = da_root, exp_out = exp_out)
    else:
        print(f"[{exp['id']}] cfg.source.enabled = false - skipping.")


def apply_features_to_source(src: dict, feat: dict) -> dict:
    """Map the shared features onto the source.* keys read by source_core."""
    _d(src, "whole", {})
    _set_from_features(src["whole"], "tmin", feat["windows"]["whole_sec"][0], "source.whole.tmin")
    _set_from_features(src["whole"], "tmax", feat["windows"]["whole_sec"][1], "source.whole.tmax")
    _d(src, "prestim", {})
    _set_from_features(src["prestim"], "tmin", feat["windows"]["pre_sec"][0], "source.prestim.tmin")
    _set_from_features(src["prestim"], "tmax", feat["windows"]["pre_sec"][1], "source.prestim.tmax")

    _d(src, "poststim", {})
    _set_from_features(src["poststim"], "tmin", feat["windows"]["post_sec"][0], "source.poststim.tmin")
    _set_from_features(src["poststim"], "tmax", feat["windows"]["post_sec"][1], "source.poststim.tmax")
    _set_from_features(src["poststim"], "phase_ref_t", feat["phase"]["post_ref_sec"], "source.poststim.phase_ref_t")

    _d(src, "lep", {})
    _set_from_features(src["lep"], "n2_window", feat["lep"]["n2_window_sec"], "source.lep.n2_window")
    _set_from_features(src["lep"], "p2_window", feat["lep"]["p2_window_sec"], "source.lep.p2_window")

    _d(src, "spectral", {})
    sp = src["spectral"]
    sp.pop("quiet_band", None)
    _set_from_features(sp, "alpha_band",         feat["bands"]["alpha_hz"],          "source.spectral.alpha_band")
    _set_from_features(sp, "slow_alpha_band",    feat["bands"]["slow_hz"],           "source.spectral.slow_alpha_band")
    _set_from_features(sp, "fast_alpha_band",    feat["bands"]["fast_hz"],           "source.spectral.fast_alpha_band")
    _set_from_features(sp, "fmin",               feat["psd"]["fmin_hz"],             "source.spectral.fmin")
    _set_from_features(sp, "fmax",               feat["psd"]["fmax_hz"],             "source.spectral.fmax")
    _set_from_features(sp, "psd_window_sec",     feat["psd"]["window_sec"],          "source.spectral.psd_window_sec")
    _set_from_features(sp, "psd_overlap_frac",   feat["psd"]["overlap_frac"],        "source.spectral.psd_overlap_frac")
    _set_from_features(sp, "psd_df_target_hz",   feat["psd"]["df_target_hz"],        "source.spectral.psd_df_target_hz")
    _set_from_features(sp, "filter_trans_bw_hz", feat["band_power"]["trans_bw_hz"],  "source.spectral.filter_trans_bw_hz")
    _set_from_features(sp, "eps0_frac",          feat["erd"]["eps0_frac"],           "source.spectral.eps0_frac")
    _set_from_features(sp, "p5_percentile",      feat["erd"]["p5_percentile"],       "source.spectral.p5_percentile")

    src["fooof"] = {k: v for k, v in feat["fooof"].items()}
    return src
