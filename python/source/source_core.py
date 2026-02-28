from __future__ import annotations

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import mne

from fooof import FOOOF

# ------------
# IO Helpers
# ------------
def src_find_preproc_set(da_root: str, exp_out: str, stage_dir: str, out_prefix: str, sub: int, allow_fallback: bool) -> str:
    """
    Pattern:
        /preproc/08_base/<out_prefix>_<sub-XXX>_fir_notch60_rs500_reref_initrej_ica_iclabel_epoch_base.set

    Root:
        <da_root>/<exp_out>/preproc/<stage_dir>
    
    We search the stage_dir for either:
        <stage_dir>/<out_prefix>*sub-XXX*_base.set
    and require exactly one match unless fallback is enabled.
    """
    sub_str = f"sub-{sub:03d}"
    preproc_root = os.path.join(da_root, exp_out, "preproc", stage_dir)

    # strict glob
    pat = os.path.join(preproc_root, f"{out_prefix}*{sub_str}*_base.set")
    hits = glob.glob(pat)

    if len(hits) == 1:
        return hits[0]
    
    if len(hits) == 0 and allow_fallback:
        # fallback: any .set containing sub-XXX and ending with _base.set
        pat2 = os.path.join(preproc_root, f"*{sub_str}*_base.set")
        hits = glob.glob(pat2)

    if len(hits) == 0:
        raise FileNotFoundError(f"No .set found for {sub_str} in {preproc_root} using pattern: {pat}")

    # choose newest
    hits.sort(key = lambda p: os.path.getmtime(p), reverse = True)
    return hits[0]

def src_read_epochs_eeglab_set(set_path: str) -> mne.Epochs:
    if hasattr(mne.io, "read_epochs_eeglab"):
        ep = mne.io.read_epochs_eeglab(set_path, verbose = "ERROR")
        ep.load_data()
        return ep
    if hasattr(mne, "read_epochs_eeglab"):
        ep = mne.read_epochs_eeglab(set_path, verbose = "ERROR")
        ep.load_data()
        return ep
    raise RuntimeError("MNE lacks read_epochs_eeglab(). Upgrade MNE or export epochs to FIF.")

# -------------
# ROI Helpers
# -------------
def src_load_labels(subjects_dir: str, parc: str):
    labels = mne.read_labels_from_annot(
        subjects = "fsaverage", parc = parc, subjects_dir = subjects_dir, verbose = "ERROR"
    )
    by_name = {lab.name: lab for lab in labels}
    return labels, by_name

def src_build_custom_rois(by_name: dict, custom_rois: dict):
    """
    custom_rois example:
        "lTC": ["superiortempotal-lh",...] (note MNE/DK uses "-lh"/"-rh" endings)
    You can supply these either:
        - MNE style: "superiortemporal-lh"
        - My style: "lh-superiortemporal"
    We'll normalize.
    """
    def normalize(lbl: str) -> str:
        s = lbl.strip()
        if s.startswith("lh-"):
            return s.replace("lh-", "") + "-lh"
        if s.startswith("rh-"):
            return s.replace("rh-", "") + "-rh"
        return s
    
    macro_labels = []
    macro_names = []

    for macro, parts in custom_rois.items():
        parts_n = [normalize(p) for p in parts]
        labs = []
        missing = []
        for p in parts_n:
            if p in by_name:
                labs.append(by_name[p])
            else:
                missing.append(p)
        if missing:
            print(f"[WARN] ROI {macro}: missing labels: {missing}")
        if not labs:
            print(f"[WARN] ROI {macro}: no labels found, skipping.")
            continue

        combined = labs[0]
        for lab in labs[1:]:
            combined = combined + lab
        combined.name = macro

        macro_labels.append(combined)
        macro_names.append(macro)

    return macro_labels, macro_names

# --------------------------
# Spectral + alpha metrics
# --------------------------
def src_psd_welch_1d(x: np.ndarray, sfreq: float, fmin: float, fmax: float):
    psd, freqs = mne.time_frequency.psd_array_welch(
        x, sfreq = sfreq, fmin = fmin, fmax = fmax, average = "mean", verbose = "ERROR"
    )
    return freqs, psd

def src_bandpower(freqs: np.ndarray, psd: np.ndarray, lo: float, hi: float) -> float:
    idx = (freqs >= lo) & (freqs <= hi)
    if not np.any(idx):
        return float("nan")
    return float(np.trapezoid(psd[idx], freqs[idx]))

def src_cog(freqs: np.ndarray, psd: np.ndarray, lo: float, hi: float) -> float:
    idx = (freqs >= lo) & (freqs <= hi)
    f = freqs[idx]
    p = psd[idx]
    denom = np.sum(p)
    if denom <= 0 or len(f) == 0:
        return float("nan")
    return float(np.sum(f * p) / denom)

def src_alpha_metrics(freqs, psd, slow, fast, alpha):
    slow_p = src_bandpower(freqs, psd, slow[0], slow[1])
    fast_p = src_bandpower(freqs, psd, fast[0], fast[1])
    alpha_p = src_bandpower(freqs, psd, alpha[0], alpha[1])
    paf = src_cog(freqs, psd, alpha[0], alpha[1])

    eps = 1e-20
    out = {
        "slow_alpha_power": slow_p,
        "fast_alpha_power": fast_p,
        "alpha_power": alpha_p,
        "paf_cog_hz": paf,
        "sf_ratio": (slow_p / fast_p) if (fast_p not in [0, np.nan] and fast_p > 0) else np.nan,
        "sf_log_ratio": float(np.log(slow_p + eps) / np.log(fast_p + eps)),
        "sf_diff": float(slow_p - fast_p),
        "sf_normdiff": float((slow_p - fast_p) / (slow_p + fast_p + eps)),
        "sf_frac_slow": float(slow_p / (slow_p + fast_p + eps)),
        "sf_frac_fast": float(fast_p / (slow_p + fast_p + eps)),
    }
    return out

# -------------
# FOOOF
# -------------
def src_run_fooof(freqs: np.ndarray, psd: np.ndarray, fooof_cfg: dict):
    fm = FOOOF(
        aperiodic_mode = fooof_cfg.get("aperiodic_mode", "fixed"),
        peak_width_limits = tuple(fooof_cfg.get("peak_width_limits", [1.0, 12.0])),
        max_n_peaks = int(fooof_cfg.get("man_n_peaks", 6)),
        min_peak_height = float(fooof_cfg.get("min_peak_height", 0.1)),
        peak_threshold = float(fooof_cfg.get("peak_threshold", 2.0)),
        verbose = False,
    )
    fr = fooof_cfg.get("freq_range", [1.0, 40.0])
    fm.fit(freqs, psd, fr)

    ap = fm.get_params("aperiodic_params") # [offset, exponent] or [offset, knee, exponent]
    out = {}

    if fooof_cfg.get("aperiodic_mode", "fixed") == "fixed":
        out["fooof_offset"] = float(ap[0])
        out["fooof_exponent"] = float(ap[1])
        out["fooof_knee"] = np.nan
    else:
        out["fooof_offset"] = float(ap[0])
        out["fooof_knee"] = float(ap[1])
        out["fooof_exponent"] = float(ap[2])

    # alpha peak extraction (best effort)
    # We'll pick the strongest peak whose CF is inside 8-12
    peaks = fm.get_params("peak_params")
    alpha_peaks = []
    for (cf, pw, bw) in peaks:
        if 8.0 <= cf <= 12.0:
            alpha_peaks.append((cf, pw, bw))
    if alpha_peaks:
        # strongest by pw
        alpha_peaks.sort(key = lambda t: t[1], reverse = True)
        cf, pw, bw = alpha_peaks[0]
        out["fooof_alpha_cf"] = float(cf)
        out["fooof_alpha_pw"] = float(pw)
        out["fooof_alpha_bw"] = float(bw)
    else:
        out["fooof_alpha_cf"] = np.nan
        out["fooof_alpha_pw"] = np.nan
        out["fooof_alpha_bw"] = np.nan

    return out, fm

def src_save_fooof_plot(path_png: str, fm: FOOOF, title: str):
    fig = fm.plot(plot_peaks = "shade", add_legend = True)
    fig.axes[0].set_title(title)
    fig.tight_layout()
    fig.savefig(path_png, dpi = 150)
    plt.close(fig)

# ----------------
# Figure helpers
# ----------------
def src_corr_matrix(X: np.ndarray) -> np.ndarray:
    """
    X: (n_rois, n_times)
    returns: (n_rois, n_rois) Pearson correlation
    """
    if X.ndim != 2:
        raise ValueError("X must be 2D (n_rois, n_times)")
    return np.corrcoef(X)

def src_save_corr_csv(path_csv: str, C: np.ndarray, roi_names: list[str]):
    df = pd.DataFrame(C, index = roi_names, columns = roi_names)
    df.to_csv(path_csv)

def src_save_corr_png(path_png: str, C: np.ndarray, roi_names: list[str], title: str):
    plt.figure(figsize = (8, 7))
    im = plt.imshow(C, vmin = 1, vmax = 1, interpolation = "nearest", aspect = "equal")
    plt.colorbar(im, fraction = 0.046, psd = 0.04)
    plt.xticks(range(len(roi_names)), roi_names, rotation = 90, fontsize = 8)
    plt.yticks(range(len(roi_names)), roi_names, fontsize = 8)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path_png, dpi = 150)
    plt.close()

def src_save_stc_snapshot_matplotlib(stc, subjects_dir: str, out_png: str, t_sec: float, title: str):
    """
    Saves a 2D snapshot using MNE's matplotlib backend.
    """
    brain = stc.plot(
        subject = "fsaverage",
        subjects_dir = subjects_dir,
        initial_time = t_sec,
        hemi = "split",
        views = "lat",
        backend = "matplotlib",
        time_viewer = False,
        colorbar = True,
        size = (900, 400),
        show = False,
    )
    # brain is a matplotlib figure for backend = "matplotlib"
    brain.suptitle(title)
    brain.savefig(out_png, dpi = 150)
    plt.close(brain)

# -----------------
# Core Execution
# -----------------
def source_core(cfg: dict, da_root: str, exp_out: str):
    src_cfg = cfg["source"]
    exp_cfg = cfg["exp"]

    out_prefix = str(exp_cfg.get("out_prefix", ""))
    if not out_prefix:
        raise ValueError("cfg.exp.out_prefix is required for strict filename resolving.")
    
    stage_dir = src_cfg["input"]["stage_dir"]
    allow_fallback = bool(src_cfg["input"]["allow_fallback_search"])
    out_root = src_cfg["outputs"]["root"]
    os.makedirs(out_root, exist_ok = True)
    
    subjects_dir = src_cfg["fsaverage"]["subjects_dir"]
    os.environ["SUBJECTS_DIR"] = subjects_dir

    parc = src_cfg["roi"]["parcellation"]
    roi_mode = "custom" if src_cfg["roi"].get("use_custom_rois", True) else "aparc_all"
    roi_extract_mode = src_cfg["roi"].get("mode", "mean_flip")

    # spectral params
    fmin =  float(src_cfg["spectra"]["fmin"])
    fmax =  float(src_cfg["spectra"]["fmax"])
    slow = tuple(src_cfg["spectra"]["slow_alpha_band"])
    fast = tuple(src_cfg["spectra"]["fast_alpha_band"])
    alpha = tuple(src_cfg["spectra"]["alpha_band"])

    # inverse params
    method = src_cfg["inverse"]["method"]
    if method.lower() != "sloreta":
        raise ValueError(f"Unsupported inverse method: {method}. Only sLORETA is currently supported.")
    snr = float(src_cfg["inverse"].get("snr", 3.0))
    lambda2 = 1.0 / (snr ** 2)

    pick_ori = src_cfg["inverse"].get("pick_ori", None)
    if isinstance(pick_ori, str) and pick_ori.lower() == "none":
        pick_ori = None

    # Noise cov window
    tmin = float(src_cfg["noise_cov"]["tmin"])
    tmax = float(src_cfg["noise_cov"]["tmax"])

    # fsaverage BEM / trans / src
    # we assume these exist under subjects_dir/fsaverage/bem
    bem_sol = os.path.join(subjects_dir, "fsaverage", "bem", "fsaverage-5120-5120-5120-bem-sol.fif")
    trans = os.path.join(subjects_dir, "fsaverage", "bem", "fsaverage-trans.fif")
    src_space = os.path.join(subjects_dir, "fsaverage", "bem", "fsaverage-ico-5-src.fif")
    
    for p in [bem_sol, trans, src_space]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing fsaverage asset: {p}")
        
    # labels
    labels_all, by_name = src_load_labels(subjects_dir, parc)

    if roi_mode == "custom":
        labels, roi_names = src_build_custom_rois(by_name, src_cfg["roi"]["custom_rois"])
    else:
        labels = labels_all
        roi_names = [l.name for l in labels]

    fooof_cfg = src_cfg.get("fooof", {})
    do_fooof = bool(fooof_cfg.get("enabled", True))
    fooof_mode = str(fooof_cfg.get("mode", "ga_only")).lower()

    for sub in exp_cfg["subjects"]:
        sub = int(sub)
        sub_str = f"sub-{sub:03d}"
        print(f"\n===== SOURCE START {sub_str} =====")

        sub_out = os.path.join(out_root, sub_str)
        fig_dir = os.path.join(sub_out, "figures")
        os.makedirs(fig_dir, exist_ok = True)

        set_path = src_find_preproc_set(da_root, exp_out, stage_dir, out_prefix, sub, allow_fallback)
        print(f"[INPUT] {set_path}")

        epochs = src_read_epochs_eeglab_set(set_path)
        sfreq = float(epochs.info["sfreq"])
        print(f"[EPOCHS] n = {len(epochs)} sfreq = {sfreq:.2f} t = ({epochs.tmin:.3f}, {epochs.tmax:.3f})")

        noise_cov = mne.compute_covariance(
            epochs, tmin = tmin, tmax = tmax, method = "empirical", rank = None, verbose = "ERROR"
        )

        fwd = mne.make_forward_solution(
            info = epochs.info,
            trans = trans,
            src = src_space,
            bem = bem_sol,
            eeg = True,
            meg = False,
            mindist = float(src_cfg["forward"].get("mindist_mm", 5.0)),
            verbose = "ERROR",
        )

        inv = mne.minimum_norm.make_inverse_operator(
            epochs.info, fwd, noise_cov, loose = 0.2, depth = 0.8, verbose = "ERROR"
        )

        stcs = mne.minimum_norm.apply_inverse_epochs(
            epochs,
            inv,
            lambda2 = lambda2,
            method = "sLORETA",
            pick_ori = pick_ori,
            verbose = "ERROR",
        )

        tc = mne.extract_label_time_course(
            stcs, labels = labels, src = fwd["src"], mode = roi_extract_mode, verbose = "ERROR"
        )
        tc = np.asarray(tc) # (n_epochs, n_rois, n_times)

        # ----- trial-wise features -----
        trial_rows = []
        fooof_trial_rows = []

        for ei in range(tc.shape[0]):
            for ri, roi in enumerate(roi_names):
                x = tc[ei, ri, :]
                freqs, psd = src_psd_welch_1d(x, sfreq, fmin, fmax)
                am = src_alpha_metrics(freqs, psd, slow, fast, alpha)

                row = {"subject": sub, "trial": ei + 1, "roi": roi, "sfreq": sfreq, **am}
                trial_rows.append(row)
                
                if do_fooof and fooof_mode == "trial_and_ga":
                    f_out, fm = src_run_fooof(freqs, psd, fooof_cfg)
                    fooof_trial_rows.append({"subeject": sub, "trial": ei + 1, "roi": roi, **f_out})

        trial_df = pd.DataFrame(trial_rows)
        trial_csv = os.path.join(sub_out, f"{sub_str}_trialwise_roi_source_spectral.csv")
        trial_df.to_csv(trial_csv, index = False)
        print(f"[SAVE] {trial_csv}")

        if do_fooof and fooof_mode == "trial_and_ga":
            f_trial_df = pd.DataFrame(fooof_trial_rows)
            f_trial_csv = os.path.join(sub_out, f"{sub_str}_trialwise_roi_fooof.csv")
            f_trial_df.to_csv(f_trial_csv, index = False)
            print(f"[SAVE] {f_trial_csv}")

        # ----- GA features -----
        ga_rows = []
        fooof_ga_rows = []

        psd_by_roi = {}
        for ri, roi in enumerate(roi_names):
            x_ga = np.mean(tc[:, ri, :], axis = 0)
            freqs, psd = src_psd_welch_1d(x_ga, sfreq, fmin, fmax)
            psd_by_roi[roi] = psd

            am = src_alpha_metrics(freqs, psd, slow, fast, alpha)
            ga_rows.append({"subject": sub, "roi": roi, "sfreq": sfreq, **am})

            if do_fooof:
                f_out, fm = src_run_fooof(freqs, psd, fooof_cfg)
                fooof_ga_rows.append({"subject": sub, "roi": roi, **f_out})

                # optional per ROI FOOOF plot (can be many; keep it GA only)
                if src_cfg["gc"].get("save_plots", True) and roi_mode == "custom":
                    png = os.path.join(fig_dir, f"{sub_str}_GA_{roi}_fooof.png")
                    src_save_fooof_plot(png, fm, title = f"{sub_str} GA {roi} FOOOF")

        ga_df = pd.DataFrame(ga_rows)
        ga_csv = os.path.join(sub_out, f"{sub_str}_GA_roi_source_spectral.csv")
        ga_df.to_csv(ga_csv, index = False)
        print(f"[SAVE] {ga_csv}")

        if do_fooof:
            f_ga_df = pd.DataFrame(fooof_ga_rows)
            f_ga_csv = os.path.join(sub_out, f"{sub_str}_GA_roi_fooof.csv")
            f_ga_df.to_csv(f_ga_csv, index = False)
            print(f"[SAVE] {f_ga_csv}")

        # QC: GA PSD overlay (custom ROIs only, otherwise too many lines)
        if src_cfg["qc"].get("save_plots", True) and roi_mode == "custom":
            plt.figure(figsize = (10, 6))
            for roi in roi_names:
                plt.plot(freqs, 10 * np.log10(psd_by_roi[roi] + 1e-20), label = roi)
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Power (dB)')
            plt.title(f"{sub_str} GA ROI PSD (sLORETA)")
            plt.legend(fontsize = 8, ncol = 2)
            plt.tight_layout()
            out_png = os.path.join(fig_dir, f"{sub_str}_GA_roi_psd.png")
            plt.savefig(out_png, dpi = 150)
            plt.close()
            print(f"[SAVE] {out_png}")

        # QC: ROI correlation matrix (custom ROIs only)
        save_conn = bool(src_cfg["outputs"].get("save_connectivity", True))
        if save_conn:
            conn_dir = os.path.join(sub_out, "connectivity")
            os.makedirs(conn_dir, exist_ok = True)

            # ---- Trial-wise corr ----
            for ei in range(tc.shape[0]):
                X = tc[ei, :, :] # (n_rois, n_times)
                C = src_corr_matrix(X)

                csv_path = os.path.join(conn_dir, f"{sub_str}_trial-{ei+1:04d}_roi_corr.csv")
                png_path = os.path.join(conn_dir, f"{sub_str}_trial-{ei+1:04d}_roi_corr.png")

                src_save_corr_csv(csv_path, C, roi_names)
                if src_cfg["qc"].get("save_plots", True):
                    src_save_corr_png(png_path, C, roi_names, title = f"{sub_str} Trial {ei+1} ROI Corr (sLORETA)")

            # ---- GA corr (mean over trials) ----
            Xga = np.mean(tc, axis = 0) # (n_rois, n_times)
            Cga = src_corr_matrix(Xga)

            csv_path = os.path.join(conn_dir, f"{sub_str}_GA_roi_corr.csv")
            png_path = os.path.join(conn_dir, f"{sub_str}_GA_roi_corr.png")

            src_save_corr_csv(csv_path, Cga, roi_names)
            if src_cfg["qc"].get("save_plots", True):
                src_save_corr_png(png_path, Cga, roi_names, title = f"{sub_str} GA ROI Corr (sLORETA)")

        save_brain = bool(src_cfg["qc"].get("save_brain_images", False)) or bool(src_cfg["outputs"].get("save_brain_images", False))
        t_snap = float(src_cfg["qc"].get("brain_snapshot_time_sec", 0.05))

        if save_brain:
            # GA/evoked STC
            evoked = epochs.average()
            stc_ga = mne.minimum_norm.apply_inverse(
                evoked, inv, lambda2 = lambda2, method = "sLORETA", pick_ori = pick_ori, verbose = "ERROR"
            )

            img_dir = os.path.join(sub_out, "brain_images")
            os.makedirs(img_dir, exist_ok = True)
            out_png = os.path.join(img_dir, f"{sub_str}_GA_sLORETA_t{t_snap:.3f}s_png")

            try:
                src_save_stc_snapshot_matplotlib(
                    stc_ga,
                    subjects_dir = subjects_dir,
                    out_png = out_png,
                    t_sec = t_snap,
                    title = f"{sub_str} GA sLORETA @ {t_snap:.3f}s",
                )
                print(f"[SAVE] {out_png}")
            except Exception as e:
                print(f"[WARN] Brain snapshot failed (matplotlib): {e}")

        print(f"===== SOURCE DONE {sub_str} =====")