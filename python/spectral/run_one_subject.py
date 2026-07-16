"""
run_one_subject.py  –  Process exactly ONE subject's signed-orientation
sLORETA ERP (pick_ori="normal"), restricted to 4 target ROIs, save checkpoint.
V 1.0.0

Designed to be invoked as a fresh subprocess per subject (see
driver_run_all_subjects.py) so the OS fully reclaims memory between subjects,
avoiding the OOM "Killed" issue seen when looping over all subjects within a
single long-running process.

Usage
-----
    python run_one_subject.py <exp_name> <sub_int>

Output
------
    <checkpoint_dir>/<exp_name>_sub-XXX_n2p2_signed.npz
        contains: ga_tc (n_rois, n_times), times (n_times,)

Skips and exits cleanly (no error) if the checkpoint already exists, so the
driver can be re-run safely after a crash without redoing finished subjects.
"""

from __future__ import annotations

import os
import sys
import glob
import traceback

import numpy as np
import mne


# =============================================================================
# CONFIG — must match gen_n2p2_erp_signed.py / exp01.json
# =============================================================================

DA_ROOT = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis"
SUBJECTS_DIR = "/home/UWO/darsenea/mne_data/"

ROIS_OF_INTEREST = ["postcentral", "supramarginal"]
TARGET_LABEL_NAMES = [f"{roi}-{hemi}" for roi in ROIS_OF_INTEREST for hemi in ("lh", "rh")]

PARCELLATION = "aparc"
ROI_MODE = "mean_flip"

SNR = 3.0
LAMBDA2 = 1.0 / (SNR ** 2)
LOOSE = 0.2
DEPTH = 0.8
PICK_ORI = "normal"

NOISE_TMIN = -0.2
NOISE_TMAX = 0.0

CHECKPOINT_DIR = os.path.join(DA_ROOT, "n2p2_signed_checkpoints")

# =============================================================================
# HELPERS
# =============================================================================

def _find_epochs_set(da_root: str, exp_name: str, sub: int) -> str | None:
    sub_str = f"sub-{sub:03d}"
    pattern = os.path.join(da_root, exp_name, "preproc", sub_str, "08_base", "*.set")
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None


def _find_fwd(da_root: str, exp_name: str, sub: int) -> str | None:
    sub_str = f"sub-{sub:03d}"
    path = os.path.join(da_root, exp_name, "source", sub_str, "fwd", f"{sub_str}_fwd.fif")
    return path if os.path.exists(path) else None


# =============================================================================
# MAIN
# =============================================================================

def main():
    if len(sys.argv) != 3:
        print("Usage: python run_one_subject.py <exp_name> <sub_int>")
        sys.exit(2)

    exp_name = sys.argv[1]
    sub = int(sys.argv[2])
    sub_str = f"sub-{sub:03d}"

    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    out_path = os.path.join(CHECKPOINT_DIR, f"{exp_name}_{sub_str}_n2p2_signed.npz")

    if os.path.exists(out_path):
        print(f"[SKIP] {exp_name}/{sub_str} — checkpoint already exists")
        sys.exit(0)

    os.environ["SUBJECTS_DIR"] = SUBJECTS_DIR

    try:
        set_path = _find_epochs_set(DA_ROOT, exp_name, sub)
        fwd_path = _find_fwd(DA_ROOT, exp_name, sub)
        if set_path is None or fwd_path is None:
            print(f"[SKIP] {exp_name}/{sub_str} — missing epochs or fwd")
            sys.exit(0)

        # ── Load labels (only the 4 needed) ───────────────────────────────────
        all_labels = mne.read_labels_from_annot(
            subject="fsaverage", parc=PARCELLATION, subjects_dir=SUBJECTS_DIR, verbose="ERROR"
        )
        by_name = {lab.name: lab for lab in all_labels}
        missing = [n for n in TARGET_LABEL_NAMES if n not in by_name]
        if missing:
            raise RuntimeError(f"Target labels not found: {missing}")
        target_labels = [by_name[n] for n in TARGET_LABEL_NAMES]

        # ── Load epochs + fwd ──────────────────────────────────────────────────
        epochs = mne.io.read_epochs_eeglab(set_path, verbose="ERROR")
        epochs.set_eeg_reference("average", projection=True, verbose="ERROR")
        fwd = mne.read_forward_solution(fwd_path, verbose="ERROR")

        noise_cov = mne.compute_covariance(
            epochs, tmin=NOISE_TMIN, tmax=NOISE_TMAX,
            method="empirical", rank=None, verbose="ERROR",
        )

        inv = mne.minimum_norm.make_inverse_operator(
            epochs.info, fwd, noise_cov,
            loose=LOOSE, depth=DEPTH, verbose="ERROR",
        )

        # ── Stream trial-by-trial, running mean ───────────────────────────────
        stcs_gen = mne.minimum_norm.apply_inverse_epochs(
            epochs, inv, lambda2=LAMBDA2, method="sLORETA",
            pick_ori=PICK_ORI, return_generator=True, verbose="ERROR",
        )

        running_sum = None
        n_trials = 0
        times = epochs.times

        for stc in stcs_gen:
            trial_tc = mne.extract_label_time_course(
                stc, labels=target_labels, src=fwd["src"],
                mode=ROI_MODE, verbose="ERROR",
            )
            trial_tc = np.asarray(trial_tc)

            if running_sum is None:
                running_sum = np.zeros_like(trial_tc)
            running_sum += trial_tc
            n_trials += 1

        if n_trials == 0:
            print(f"[SKIP] {exp_name}/{sub_str} — no trials processed")
            sys.exit(0)

        ga_tc = running_sum / n_trials   # (n_rois, n_times)

        np.savez(out_path, ga_tc=ga_tc, times=times, n_trials=n_trials)
        print(f"[OK] {exp_name}/{sub_str} -> {out_path}  (n_trials={n_trials})")
        sys.exit(0)

    except Exception:
        print(f"[ERROR] {exp_name}/{sub_str}:\n{traceback.format_exc()}")
        sys.exit(1)


if __name__ == "__main__":
    main()