"""
source_core.py  -  Per-subject orchestration loop for sLORETA source analysis.
V 4.1.0

V4.0.0 changes vs V3.1.0 (channel/source parity release; MATLAB twin is
spectral_core.m V3.0.0):
    - All feature estimators come from the shared definitions in
      src_spectral.py V3.0.0 (explicit FIR + Hilbert band power, explicit
      Welch, shared window masks), which have line-for-line MATLAB twins.
    - All feature parameters come from cfg["features"] via source_default.py.
    - ERD guard eps0 is data-relative (no 45-55 Hz quiet band).
    - Phase / ITC use the same FIR as band power.
    - LEP searched on the full epoch with shared masks; windows from features.
    - FOOOF on the subject x ROI trial-mean whole-epoch PSD, alpha band from
      features, now also writes fooof_r2 / fooof_error.
    - roi.hemisphere_split adds <ROI>_lh / <ROI>_rh next to each bilateral ROI.
    - V4.1.0: whole window from features.windows.whole_sec (not the full
      epoch); slow_alpha_frac dropped; delta_erd -> erd_asym.

Per-subject outputs
───────────────────
    <out_root>/sub-XXX/csv/
        sub-XXX_source_trial.csv      trial x ROI      (read by the R merge)
        sub-XXX_source_ga.csv         subject x ROI
        sub-XXX_source_ga_fooof.csv   subject x ROI FOOOF only
        sub-XXX_source_ga_timecourse.npy, sub-XXX_source_times.npy
    <out_root>/sub-XXX/fwd/sub-XXX_fwd.fif   cached forward solution
    <out_root>/sub-XXX/logs/

Column names in sub-XXX_source_trial.csv are identical to the metric columns
of MATLAB's sub-XXX_spectral_chan_by_trial.csv.
"""

from __future__ import annotations

import os
import traceback

import numpy as np
import pandas as pd

from src_io       import src_open_log, src_logmsg, src_close_log, src_find_set, src_read_epochs
from src_assets   import src_load_fsaverage_assets, src_load_labels, src_build_custom_rois
from src_inverse  import src_make_inverse_operator, src_apply_inverse_epochs
from src_alpha_features import src_compute_band_power_tcs, src_compute_window_alpha_features
from src_add_prefix import src_add_prefix_rows
from src_merge_rows import src_merge_rows
from src_compute_metric_deltas import src_compute_metric_deltas
from src_erd      import src_compute_noise_stats, src_compute_erd_metrics
from src_prestim  import src_compute_prestim_phase, src_compute_ga_prestim_phase, src_compute_tvi_alpha
from src_poststim import src_compute_poststim_phase, src_compute_ga_poststim_phase, src_compute_itc
from src_lep      import src_compute_lep_trial, src_compute_lep_ga
from src_fooof    import fooof_available, fooof_package_name, src_compute_fooof_ga
from src_write    import src_write_trial_csv, src_write_ga_csv, src_write_fooof_csv
from src_plot     import src_plot_ga_timecourse


def compute_subject_features(tc: np.ndarray, times: np.ndarray, sfreq: float,
                             src_cfg: dict, n_rois: int, sub: int, logf = None):
    """
    All source features for one subject from ROI time courses.

    Kept separate from I/O and the inverse so tests/parity_test.py can run
    it on synthetic data.

    Args:
        tc     : (n_epochs, n_rois, n_times) ROI time courses, full epoch
        times  : (n_times,) seconds
        sfreq  : Hz
        src_cfg: cfg["source"] after source_default()

    Returns:
        trial_rows, ga_rows, fooof_rows  (lists of dicts keyed by roi_idx)
    """
    log = (lambda *a: src_logmsg(logf, *a)) if logf is not None else (lambda *a: None)
    sp = src_cfg["spectral"]

    fmin, fmax = float(sp["fmin"]), float(sp["fmax"])
    alpha = tuple(sp["alpha_band"])
    slow = tuple(sp["slow_alpha_band"])
    fast = tuple(sp["fast_alpha_band"])
    trans_bw = float(sp["filter_trans_bw_hz"])
    psd_kw = dict(psd_window_sec = float(sp["psd_window_sec"]),
                  df_target = float(sp["psd_df_target_hz"]),
                  overlap = float(sp.get("psd_overlap_frac", 0.5)))
    eps0_frac = float(sp.get("eps0_frac", 1e-3))
    p_pct = float(sp.get("p5_percentile", 5.0))

    pre_tmin, pre_tmax = float(src_cfg["prestim"]["tmin"]), float(src_cfg["prestim"]["tmax"])
    post_tmin, post_tmax = float(src_cfg["poststim"]["tmin"]), float(src_cfg["poststim"]["tmax"])
    post_ref_t = float(src_cfg["poststim"]["phase_ref_t"])
    n2_window = tuple(src_cfg["lep"]["n2_window"])
    p2_window = tuple(src_cfg["lep"]["p2_window"])
    fooof_cfg = src_cfg.get("fooof", {})

    # Band-power time courses (full epoch) + whole / pre / post features
    log("[BANDPOW] Filter-Hilbert power: slow %s, fast %s, alpha %s Hz (trans_bw = %.2f Hz)",
        str(slow), str(fast), str(alpha), trans_bw)
    power_tcs = src_compute_band_power_tcs(tc, sfreq, alpha, slow, fast, trans_bw)

    win = lambda a, b: src_compute_window_alpha_features(power_tcs, tc, times, a, b, sfreq, alpha, fmin, fmax, **psd_kw)
    whole_rows, ga_whole_rows, psd_by_roi_whole = win(float(src_cfg["whole"]["tmin"]), float(src_cfg["whole"]["tmax"]))
    pre_rows,   ga_pre_rows,   _ = win(pre_tmin, pre_tmax)
    post_rows,  ga_post_rows,  _ = win(post_tmin, post_tmax)

    delta_rows    = src_compute_metric_deltas(pre_rows, post_rows, key_cols = ("trial", "roi_idx"))
    ga_delta_rows = src_compute_metric_deltas(ga_pre_rows, ga_post_rows, key_cols = ("roi_idx",))

    # ERD + p5_flag (per-ROI eps0 / thresholds)
    noise_stats = src_compute_noise_stats(pre_rows, eps0_frac, p_pct)
    erd_rows    = src_compute_erd_metrics(pre_rows, post_rows, noise_stats, ("trial", "roi_idx"), use_p5_flag = True)
    ga_erd_rows = src_compute_erd_metrics(ga_pre_rows, ga_post_rows, noise_stats, ("roi_idx",), use_p5_flag = False)

    # Phase (full-epoch filtered signal) + ITC
    phase_pre_rows  = src_compute_prestim_phase(tc, times, sfreq, slow, trans_bw)
    phase_post_rows = src_compute_poststim_phase(tc, times, sfreq, slow, post_ref_t, trans_bw)
    tc_ga = np.mean(tc, axis = 0, keepdims = True)
    phase_pre_ga_rows  = src_compute_ga_prestim_phase(tc_ga, times, sfreq, slow, trans_bw)
    phase_post_ga_rows = src_compute_ga_poststim_phase(tc_ga, times, sfreq, slow, post_ref_t, trans_bw)
    itc_ga_rows = src_compute_itc(tc, times, sfreq, slow, post_tmin, post_tmax, trans_bw)

    # TVI from the per-trial pre_sf_balance sequence of each ROI
    pre_df = pd.DataFrame(pre_rows)
    tvi_rows = [{"roi_idx": ri,
                 "TVI_alpha": src_compute_tvi_alpha(pre_df.loc[pre_df["roi_idx"] == ri].sort_values("trial")["sf_balance"].values)}
                for ri in range(n_rois)]

    # LEP
    lep_trial_rows = src_compute_lep_trial(tc, times, sfreq, n2_window, p2_window)
    ga_lep_rows    = src_compute_lep_ga(tc, times, sfreq, n2_window, p2_window)

    # FOOOF on the subject x ROI trial-mean whole-epoch PSD
    fooof_rows: list = []
    if bool(fooof_cfg.get("enabled", True)):
        if not fooof_available():
            log("[WARN] FOOOF enabled but neither 'specparam' nor 'fooof' is installed; skipping.")
        else:
            log("[FOOOF] Fitting subject x ROI trial-mean PSDs (%s)...", fooof_package_name())
            fooof_rows, _ = src_compute_fooof_ga(psd_by_roi_whole, n_rois, sub, fooof_cfg, alpha)

    trial_rows = src_merge_rows(
        src_add_prefix_rows(whole_rows, "whole_"),
        src_add_prefix_rows(pre_rows, "pre_"),
        src_add_prefix_rows(post_rows, "post_"),
        delta_rows, erd_rows, phase_pre_rows, phase_post_rows, lep_trial_rows,
        key_cols = ("trial", "roi_idx"),
    )
    ga_rows = src_merge_rows(
        src_add_prefix_rows(ga_whole_rows, "whole_"),
        src_add_prefix_rows(ga_pre_rows, "pre_"),
        src_add_prefix_rows(ga_post_rows, "post_"),
        ga_delta_rows, ga_erd_rows, phase_pre_ga_rows, phase_post_ga_rows,
        ga_lep_rows, tvi_rows, itc_ga_rows, fooof_rows,
        key_cols = ("roi_idx",),
    )
    return trial_rows, ga_rows, fooof_rows


def source_core(cfg: dict, da_root: str, exp_out: str):
    """Per-subject sLORETA loop. Called by source_default() after validation."""
    src_cfg = cfg["source"]
    exp_cfg = cfg["exp"]

    out_prefix   = str(exp_cfg.get("out_prefix", ""))
    stage_dir    = src_cfg["input"]["stage_dir"]
    allow_fb     = bool(src_cfg["input"]["allow_fallback_search"])
    out_root     = src_cfg["outputs"]["root"]
    subjects_dir = src_cfg["fsaverage"]["subjects_dir"]
    os.environ["SUBJECTS_DIR"] = subjects_dir

    roi_cfg    = src_cfg["roi"]
    parc       = roi_cfg["parcellation"]
    use_custom = bool(roi_cfg.get("use_custom_rois", True))
    hemi_split = bool(roi_cfg.get("hemisphere_split", True))
    roi_mode   = roi_cfg.get("mode", "mean_flip")

    noise_tmin = float(src_cfg["noise_cov"]["tmin"])
    noise_tmax = float(src_cfg["noise_cov"]["tmax"])
    mindist_mm = float(src_cfg["forward"]["mindist_mm"])
    loose      = float(src_cfg["inverse"].get("loose", 0.2))
    depth      = float(src_cfg["inverse"].get("depth", 0.8))
    lambda2    = 1.0 / (float(src_cfg["inverse"].get("snr", 3.0)) ** 2)
    pick_ori_raw = src_cfg["inverse"].get("pick_ori", "normal")
    pick_ori   = None if (pick_ori_raw is None or str(pick_ori_raw).lower() == "none") else str(pick_ori_raw)

    n2_window = tuple(src_cfg["lep"]["n2_window"])
    p2_window = tuple(src_cfg["lep"]["p2_window"])
    do_brain  = bool(src_cfg["qc"].get("save_brain_images", False))

    os.makedirs(out_root, exist_ok = True)

    # ── Shared assets ────────────────────────────────────────────────────────
    bem_sol, trans, src_space = src_load_fsaverage_assets(subjects_dir)
    labels_all, by_name = src_load_labels(subjects_dir, parc)

    if use_custom and roi_cfg.get("custom_rois"):
        custom_rois = {k: v for k, v in roi_cfg["custom_rois"].items() if isinstance(v, list)}
        labels, roi_names = src_build_custom_rois(by_name, custom_rois, hemisphere_split = hemi_split)
    else:
        labels = labels_all
        roi_names = [lab.name for lab in labels]

    print(f"[LABELS] {len(labels)} ROIs: {', '.join(roi_names)}")

    # ── Subject loop ─────────────────────────────────────────────────────────
    for sub in [int(s) for s in exp_cfg["subjects"]]:
        sub_str = f"sub-{sub:03d}"
        print(f"\n{'=' * 60}\n  SOURCE START  {sub_str}\n{'=' * 60}")

        sub_out = os.path.join(out_root, sub_str)
        csv_dir = os.path.join(sub_out, "csv")
        fig_dir = os.path.join(sub_out, "figures")
        fwd_dir = os.path.join(sub_out, "fwd")
        log_dir = os.path.join(sub_out, "logs")
        for d in (sub_out, csv_dir, fig_dir, fwd_dir, log_dir):
            os.makedirs(d, exist_ok = True)

        logf = src_open_log(log_dir, sub)
        try:
            set_path = src_find_set(da_root, exp_out, stage_dir, out_prefix, sub, allow_fb)
            src_logmsg(logf, "[LOAD] %s", set_path)

            epochs = src_read_epochs(set_path)
            sfreq = float(epochs.info["sfreq"])
            src_logmsg(logf, "[EPOCHS] n=%d  sfreq=%.1f Hz  t=(%.3f, %.3f)s",
                       len(epochs), sfreq, epochs.tmin, epochs.tmax)
            if len(epochs) == 0:
                src_logmsg(logf, "[SKIP] No epochs - skipping subject.")
                continue

            fwd_cache = os.path.join(fwd_dir, f"{sub_str}_fwd.fif")
            inv, fwd = src_make_inverse_operator(
                epochs, bem_sol, trans, src_space,
                noise_tmin, noise_tmax, mindist_mm, loose, depth,
                fwd_cache_path = fwd_cache, logf = logf,
            )

            tc, times = src_apply_inverse_epochs(epochs, inv, fwd, labels, lambda2, pick_ori, roi_mode, logf)
            src_logmsg(logf, "[TC] shape: %s  (epochs x ROIs x times)", str(tc.shape))

            trial_rows, ga_rows, fooof_rows = compute_subject_features(
                tc, times, sfreq, src_cfg, len(roi_names), sub, logf)

            np.save(os.path.join(csv_dir, f"{sub_str}_source_ga_timecourse.npy"), np.mean(tc, axis = 0))
            np.save(os.path.join(csv_dir, f"{sub_str}_source_times.npy"), times)

            src_write_trial_csv(os.path.join(csv_dir, f"{sub_str}_source_trial.csv"), sub, roi_names, trial_rows, logf)
            src_write_ga_csv(os.path.join(csv_dir, f"{sub_str}_source_ga.csv"), sub, roi_names, ga_rows, logf)
            if fooof_rows:
                src_write_fooof_csv(os.path.join(csv_dir, f"{sub_str}_source_ga_fooof.csv"),
                                    sub, roi_names, fooof_rows, logf)

            if do_brain:
                src_plot_ga_timecourse(
                    os.path.join(fig_dir, f"{sub_str}_source_GA_timecourse"),
                    tc, times, roi_names, sub_str, n2_window, p2_window, logf,
                )

            src_logmsg(logf, "==== SOURCE DONE %s ====", sub_str)

        except FileNotFoundError as e:
            src_logmsg(logf, "[SKIP] %s - File not found: %s", sub_str, str(e))
        except ValueError as e:
            src_logmsg(logf, "[SKIP] %s - config/data issue: %s\n%s", sub_str, str(e), traceback.format_exc())
        except Exception:
            src_logmsg(logf, "[ERROR] %s - unexpected error:\n%s", sub_str, traceback.format_exc())
        finally:
            src_close_log(logf)
