"""
source_core.py  –  Per-subject orchestration loop for sLORETA source analysis.
V 3.0.0

V3.0.0 changes vs V2.2.0 (whole/pre/post/delta metric generalization, to
match spectral_core.m V2.2.0 -- see src_alpha_features.py, src_erd.py,
src_compute_metric_deltas.py, src_merge_rows.py, src_add_prefix.py for the
new shared modules this pulls in):
    - Every alpha feature (10 from src_compute_alpha_features, plus
    psi_cog) is now computed for THREE windows and prefixed accordingly:
    whole_<name> (full epoch), pre_<name>, post_<name>
    - Every metric additionally gets a generic delta_<name> = post - pre
    (src_compute_metric_deltas.py), on top of the ERD-specific
    functional change
    - p5_flag added
    - 4 previously-missing metric families added: rel_slow_alpha,
    rel_fast_alpha, sf_ratio, slow_alpha_frac
    - Field names now match spectral_core.m exactly (snake_case,
    whole_/pre_/post_/delta_ predixes) so channel-space and source-space
    metrics can be compared directly by metric_name

Call chain:
    exp01_source.py  ->  source_default.py  ->  source_core()

Per-subject output layout
──────────────────────────
    <out_root>/sub-XXX/
        csv/
            sub-XXX_source_trial.csv          all per-trial metrics (all ROIs)
            sub-XXX_source_ga.csv             grand-average metrics + TVI + ITC (all ROIs)
            sub-XXX_source_ga_fooof.csv       FOOOF aperiodic (if enabled, all ROIs)
        figures/
            sub-XXX_source_GA_brain_t<T>s_lh.png   GA sLORETA brain snapshot
            sub-XXX_source_GA_brain_t<T>s_rh.png
        fwd/
            sub-XXX_fwd.fif                   cached forward solution
        logs/
            sub-XXX_source_<timestamp>.log

Per-ROI figures (ERD, LEP, FOOOF, phase) are not generated at this stage.
They will be added in a later pass after R GAMMs identify significant ROIs,
so only the relevant subset is ever rendered.
"""

from __future__ import annotations

import os
import traceback

import numpy as np
import pandas as pd

from src_io       import src_open_log, src_logmsg, src_close_log, src_find_set, src_read_epochs
from src_assets   import src_load_fsaverage_assets, src_load_labels, src_build_custom_rois
from src_inverse  import src_make_inverse_operator, src_apply_inverse_epochs
from src_alpha_features import src_compute_window_alpha_features, src_compute_ga_window_alpha_features
from src_add_prefix import src_add_prefix_rows
from src_merge_rows import src_merge_rows
from src_compute_metric_deltas import src_compute_metric_deltas
from src_erd import src_compute_noise_stats, src_compute_erd_metrics
from src_prestim  import (src_compute_prestim_metrics, src_compute_ga_prestim_metrics, src_compute_tvi_alpha)
from src_poststim import (src_compute_poststim_metrics, src_compute_ga_poststim_metrics, src_compute_itc)
from src_lep      import src_compute_lep_trial, src_compute_lep_ga
from src_fooof    import fooof_available, fooof_package_name, src_compute_fooof_ga
from src_write    import src_write_trial_csv, src_write_ga_csv, src_write_fooof_csv
from src_plot     import src_plot_ga_timecourse


# =============================================================================
# HELPERS
# =============================================================================

def _crop(tc_full: np.ndarray, t_full: np.ndarray,
          tmin: float, tmax: float, label: str):
    """Crop tc_full to the [tmin, tmax] time window. Raises ValueError if empty."""
    mask = (t_full >= tmin) & (t_full <= tmax)
    if not np.any(mask):
        raise ValueError(
            f"{label} window [{tmin:.3f}, {tmax:.3f}]s has no samples "
            f"in epoch range [{t_full[0]:.3f}, {t_full[-1]:.3f}]s."
        )
    return tc_full[:, :, mask], t_full[mask]


# =============================================================================
# MAIN LOOP
# =============================================================================

def source_core(cfg: dict, da_root: str, exp_out: str):
    """
    Per-subject sLORETA source localization loop.
    Called by source_default() after config validation and path resolution.
    """
    src_cfg = cfg["source"]
    exp_cfg = cfg["exp"]

    # ── Config ────────────────────────────────────────────────────────────────
    out_prefix   = str(exp_cfg.get("out_prefix", ""))
    stage_dir    = src_cfg["input"]["stage_dir"]
    allow_fb     = bool(src_cfg["input"]["allow_fallback_search"])
    out_root     = src_cfg["outputs"]["root"]
    subjects_dir = src_cfg["fsaverage"]["subjects_dir"]
    os.environ["SUBJECTS_DIR"] = subjects_dir

    parc       = src_cfg["roi"]["parcellation"]
    use_custom = bool(src_cfg["roi"].get("use_custom_rois", False))
    roi_mode   = src_cfg["roi"].get("mode", "mean_flip")

    fmin  = float(src_cfg["spectral"]["fmin"])
    fmax  = float(src_cfg["spectral"]["fmax"])
    alpha = tuple(src_cfg["spectral"]["alpha_band"])
    slow  = tuple(src_cfg["spectral"]["slow_alpha_band"])
    fast  = tuple(src_cfg["spectral"]["fast_alpha_band"])

    pre_tmin   = float(src_cfg["prestim"]["tmin"])
    pre_tmax   = float(src_cfg["prestim"]["tmax"])
    post_tmin  = float(src_cfg["poststim"]["tmin"])
    post_tmax  = float(src_cfg["poststim"]["tmax"])
    post_ref_t = float(src_cfg["poststim"].get("phase_ref_t", 0.2))
    noise_tmin = float(src_cfg["noise_cov"]["tmin"])
    noise_tmax = float(src_cfg["noise_cov"]["tmax"])

    n2_window = tuple(src_cfg["lep"]["n2_window"])
    p2_window = tuple(src_cfg["lep"]["p2_window"])

    mindist_mm   = float(src_cfg["forward"]["mindist_mm"])
    loose        = float(src_cfg["inverse"].get("loose", 0.2))
    depth        = float(src_cfg["inverse"].get("depth", 0.8))
    snr          = float(src_cfg["inverse"].get("snr", 3.0))
    lambda2      = 1.0 / (snr ** 2)
    pick_ori_raw = src_cfg["inverse"].get("pick_ori", None)
    pick_ori     = None if (pick_ori_raw is None or
                            str(pick_ori_raw).lower() == "none") \
                       else str(pick_ori_raw)

    fooof_cfg = src_cfg.get("fooof", {})
    do_fooof  = bool(fooof_cfg.get("enabled", True))

    do_brain = bool(src_cfg["qc"].get("save_brain_images", False))

    os.makedirs(out_root, exist_ok=True)

    # ── Load shared assets once ───────────────────────────────────────────────
    bem_sol, trans, src_space = src_load_fsaverage_assets(subjects_dir)
    labels_all, by_name       = src_load_labels(subjects_dir, parc)

    if use_custom and src_cfg["roi"].get("custom_rois"):
        custom_rois = {
            k: v for k, v in src_cfg["roi"]["custom_rois"].items()
            if isinstance(v, list)
        }
        labels, roi_names = src_build_custom_rois(by_name, custom_rois)
    else:
        labels    = labels_all
        roi_names = [lab.name for lab in labels]

    print(f"[LABELS] {len(labels)} ROIs loaded "
          f"({'custom' if use_custom else parc + ' all-labels'})")

    # ── Subject loop ──────────────────────────────────────────────────────────
    for sub in [int(s) for s in exp_cfg["subjects"]]:
        sub_str = f"sub-{sub:03d}"
        print(f"\n{'='*60}\n  SOURCE START  {sub_str}\n{'='*60}")

        sub_out = os.path.join(out_root, sub_str)
        csv_dir = os.path.join(sub_out, "csv")
        fig_dir = os.path.join(sub_out, "figures")
        fwd_dir = os.path.join(sub_out, "fwd")
        log_dir = os.path.join(sub_out, "logs")
        for d in [sub_out, csv_dir, fig_dir, fwd_dir, log_dir]:
            os.makedirs(d, exist_ok=True)

        logf = src_open_log(log_dir, sub)

        try:
            # 1. Load epochs ──────────────────────────────────────────────────
            set_path = src_find_set(da_root, exp_out, stage_dir, out_prefix, sub, allow_fb)
            src_logmsg(logf, "[LOAD] %s", set_path)

            epochs = src_read_epochs(set_path)
            sfreq  = float(epochs.info["sfreq"])
            src_logmsg(logf, "[EPOCHS] n=%d  sfreq=%.1f Hz  t=(%.3f, %.3f)s",
                       len(epochs), sfreq, epochs.tmin, epochs.tmax)

            if len(epochs) == 0:
                src_logmsg(logf, "[SKIP] No epochs — skipping subject.")
                continue

            # 2. Forward solution + inverse operator ──────────────────────────
            fwd_cache = os.path.join(fwd_dir, f"{sub_str}_fwd.fif")
            inv, fwd  = src_make_inverse_operator(
                epochs, bem_sol, trans, src_space,
                noise_tmin, noise_tmax, mindist_mm, loose, depth,
                fwd_cache_path=fwd_cache, logf=logf,
            )

            # 3. Apply sLORETA to all epochs ───────────────────────────────────
            tc, times = src_apply_inverse_epochs(
                epochs, inv, fwd, labels, lambda2, pick_ori, roi_mode, logf
            )
            src_logmsg(logf, "[TC] shape: %s  (epochs × ROIs × times)", str(tc.shape))

            # 4. Crop windows ──────────────────────────────────────────────────
            tc_pre,  times_pre  = _crop(tc, times, pre_tmin,  pre_tmax,  "Pre-stim")
            tc_post, times_post = _crop(tc, times, post_tmin, post_tmax, "Post-stim")
            src_logmsg(logf, "[PRE]  %.3f – %.3f s  (%d samples)",
                       times_pre[0],  times_pre[-1],  len(times_pre))
            src_logmsg(logf, "[POST] %.3f – %.3f s  (%d samples)",
                       times_post[0], times_post[-1], len(times_post))

            # 5. Whole / pre / post alpha features ──────────────────────────────
            src_logmsg(logf, "[FEAT] Whole-epoch alpha features...")
            whole_rows = src_compute_window_alpha_features(
                tc, sfreq, alpha, slow, fast, fmin, fmax,
            )
            ga_whole_rows, psd_by_roi_whole = src_compute_ga_window_alpha_features(
                tc, sfreq, alpha, slow, fast, fmin, fmax,
            )

            src_logmsg(logf, "[FEAT] Pre-stimulus alpha features...")
            pre_rows = src_compute_window_alpha_features(
                tc_pre, sfreq, alpha, slow, fast, fmin, fmax,
            )
            ga_pre_rows, _ = src_compute_ga_window_alpha_features(
                tc_pre, sfreq, alpha, slow, fast, fmin, fmax,
            )

            src_logmsg(logf, "[FEAT] Post-stimulus alpha features...")
            post_rows = src_compute_window_alpha_features(
                tc_post, sfreq, alpha, slow, fast, fmin, fmax,
            )
            ga_post_rows, _ = src_compute_ga_window_alpha_features(
                tc_post, sfreq, alpha, slow, fast, fmin, fmax,
            )

            # Generic delta_<name> = post - pre
            delta_rows = src_compute_metric_deltas(pre_rows, post_rows, key_cols = ("trial", "roi_idx"))
            ga_delta_rows = src_compute_metric_deltas(ga_pre_rows, ga_post_rows, key_cols = ("roi_idx",))

            # ERD family (fractional change, power metrics only) + p5_flag.
            src_logmsg(logf, "[ERD] Computing ERD family (erd_slow, erd_fast, erd_pow_alpha_total, delta_erd) + p5_flag...")
            noise_stats = src_compute_noise_stats(tc_pre, sfreq, slow, fast, fmin, fmax)
            erd_rows = src_compute_erd_metrics(
                pre_rows, post_rows, noise_stats, key_cols = ("trial", "roi_idx"), use_p5_flag = True,
            )
            ga_erd_rows = src_compute_erd_metrics(
                ga_pre_rows, ga_post_rows, noise_stats, key_cols = ("roi_idx",), use_p5_flag = False,
            )

            # 6. Hilbert phase (unchanged windows: t = 0 for pre, poststim_ref_t
            #    for post; both derived from the full epoch, not the cropped
            #    pre/post windows, to avoid Hilbert edge effects) ---------------
            src_logmsg(logf, "[FEAT] Hilbert phase (pre + post)...")
            phase_pre_rows = src_compute_prestim_phase(tc, times, sfreq, slow)
            phase_post_rows = src_compute_poststim_phase(tc, times, sfreq, slow, post_ref_t)

            # 7. TVI_alpha (pre-stim-oly, unchanged in behaviour; input
            #    renamed BI_pre -> pre_sf_balance) ------------------------------
            tvi_by_roi: dict = {}
            pre_df = pd.DataFrame(pre_rows)
            for ri in range(len(roi_names)):
                bi_seq = pre_df.loc[pre_df["roi_idx"] == ri, "sf_balance"].values
                tvi_by_roi[ri] = src_compute_tvi_alpha(bi_seq)

            # 8. ITC (unchanged) + GA time courses / phase -------------------------
            itc_rows = src_compute_itc(tc, times, sfreq, slow, post_tmin, post_tmax)

            itc_ga = np.mean(tc, axis = 0, keepdims = True)

            phase_pre_ga_rows = src_compute_ga_prestim_phase(tc_ga, times, sfreq, slow)
            phase_post_ga_rows = src_compute_ga_poststim_phase(tc_ga, times, sfreq, slow, post_ref_t)

            ga_tc = np.mean(tc, axis = 0)
            npy_path = os.path.join(csv_dir, f"{sub_str}_source_ga_timecourse.npy")
            np.save(npy_path, ga_tc)
            src_logmsg(logf, "[NPY] Saved GA timecourse: %s shape = %s", npy_path, str(ga_tc.shape))

            times_path = os.path.join(csv_dir, f"{sub_str}_source_times.npy")
            np.save(times_path, times)
            src_logmsg(logf, "[NPY] Saved times: %s", times_path)

            # 9. LEP features -----------------------------------------------------
            src_logmsg(logf, "[FEAT] LEP features (N2/P2)...")
            lep_trial_rows = src_compute_lep_trial(tc_post, times_post, n2_window, p2_window)
            ga_lep_rows = src_compute_lep_ga(tc_post, times_post, n2_window, p2_window)

            # 10. FOOOF -----------------------------------------------------------
            # NOTE: Now fits the whole-epoch GA PSD, not the pre-stim GA PSD
            fooof_rows: list = []
            if do_fooof:
                if not fooof_available():
                    src_logmsg(logf,
                               "[WARN] FOOOF enabled but neither 'specparam' nor 'fooof' is installed, skipping.")
                else:
                    src_logmsg(logf, "[FOOOF] Fitting GA PSDs (%s)...", fooof_package_name())
                    fooof_rows = src_compute_fooof_ga(
                        psd_by_roi_whole, len(roi_names), sub, fooof_cfg,
                    )

            # 11. Assemble + write CSVs ------------------------------------------------
            trial_rows_merged = src_merge_rows(
                src_add_prefix_rows(whole_rows, "whole_"),
                src_add_prefix_rows(pre_rows, "pre_"),
                src_add_prefix_rows(post_rows, "post_"),
                delta_rows,
                erd_rows,
                phase_pre_rows,
                phase_post_rows,
                lep_trial_rows,
                key_cols = ("trial", "roi_idx"),
            )

            tvi_rows = [{"roi_idx": ri, "TVI_alpha": tvi_by_roi[ri]} for ri in range(len(roi_names))]

            ga_rows_merged = src_merge_rows(
                src_add_prefix_rows(ga_whole_rows, "whole_"),
                src_add_prefix_rows(ga_pre_rows, "pre_"),
                src_add_prefix_rows(ga_post_rows, "post_"),
                ga_delta_rows,
                ga_erd_rows,
                phase_pre_ga_rows,
                phase_post_ga_rows,
                ga_lep_rows,
                tvi_rows,
                itc_ga_rows,
                fooof_rows,
                key_cols = ("roi_idx",)
            )

            src_write_trial_csv(
                os.path.join(csv_dir, f"{sub_str}_source_trial.csv"),
                sub, roi_names, trial_rows_merged, logf,
            )
            src_write_ga_csv(
                os.path.join(csv_dir, f"{sub_str}_source_ga_fooof.csv"),
                sub, roi_names, ga_rows_merged, logf,
            )
            if fooof_rows:
                src_write_fooof_csv(
                    os.path.join(csv_dir, f"{sub_str}_source_ga_fooof.csv"),
                    roi_names, fooof_rows, logf,
                )

            # 12. QC heatmap (sanity check only) --------------------------------------
            if do_brain:
                src_plot_ga_timecourse(
                    os.path.join(fig_dir, f"{sub_str}_source_GA_timecourse"),
                    tc, times, roi_names, sub_str,
                    n2_window, p2_window, logf,
                )

            src_logmsg(logf, "==== SOURCE DONE %s ====", sub_str)

        except FileNotFoundError as e:
            src_logmsg(logf, "[SKIP] %s - File not found: %s", sub_str, str(e))
        except ValueError as e:
            src_logmsg(logf, "[SKIP] %s - config/data issue: %s", sub_str, str(e))
        except Exception:
            src_logmsg(logf, "[ERROR] %s - unexpected error:\n%s", sub_str, traceback.format_exc())
        finally:
            src_close_log(logf)