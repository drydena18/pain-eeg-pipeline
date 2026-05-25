"""
src_io.py - IO helpers for the source localization pipeline.

Covers:
    - Log file open / write / close
    - Preprocessed .set file resolution
    - EEGLAB .set -> MNE Epochs loading (version-compatible shim)
"""

from __future__ import annotations

import datetime
import glob
import os

import mne

# -- MNE version shim for read_epochs_eeglab -----
# The function moved between MNE versions; resolve once at import time
def _resolve_read_epochs_eeglab():
    if hasattr(mne.io, "read_epochs_eeglab"):
        return mne.io.read_epochs_eeglab
    if hasattr(mne, "read_epochs_eeglab"):
        return mne.read_epochs_eeglab
    raise RuntimeError(
        "mne.io.read_epochs_eeglab not found. "
        "Upgrade MNE (>= 1.0) or export epochs to FIF first."
    )

_read_epochs_eeglab = _resolve_read_epochs_eeglab()

# ====================================================================
# LOG HELPERS
# ====================================================================
def src_open_log(log_dir: str, sub: int, tag: str = "source"):
    """
    Open a per-subject text file (line-buffered). Returns a file handle.
    """
    os.makedirs(log_dir, exist_ok = True)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join(log_dir, f"sub-{sub:03d}_{tag}_{ts}.log")
    return open(path, "w", buffering = 1)

def src_logmsg(logf, fmt: str, *args):
    """
    Write a timestamped message to logf and stdout simultaneously.
    """
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    msg = (fmt % args) if args else fmt
    line = f"[{ts}] {msg}"
    print(line)
    if logf is not None:
        logf.write(line + "\n")

def src_close_log(logf):
    """
    Safely close a log file handle.
    """
    try:
        if logf is not None:
            logf.close()
    except Exception:
        pass

# ====================================================================
# FILE RESOLUTION
# ====================================================================
def src_find_set(
        da_root: str,
        exp_out: str,
        stage_dir: str,
        out_prefix: str,
        sub: int,
        allow_fallback: bool,
) -> str:
    """
    Locate the preprocessed .set file for a subject inside a given stage directory.

    Expected layout:
        <da_root>/<exp_out>/preproc/<sub_str>/<stage_dir>/<prefix><sub_int>*_base.set

    Example:
        .../da-analysis/26ByBiosemi/preproc/sub-001/08_base/26BB_62_001_..._base.set

    Search order:
        1. <preproc_root>/<out_prefix>*<sub_int>*_base.set  (prefix + bare int)
        2. <preproc_root>/*_base.set                         (any _base.set, if fallback)
        3. When multiple hits exist, the newest file by modification time wins.

    Args:
        da_root         : Root da-analysis directory (e.g., .../CNED/da-analysis)
        exp_out         : Experiment output folder name (e.g., "26ByBiosemi")
        stage_dir       : Preprocessing stage subdirectory (e.g., "08_base")
        out_prefix      : Filename prefix (e.g., "26BB_62_")
        sub             : Integer subject ID
        allow_fallback  : When True, retry without the out_prefix constraint

    Raises:
        FileNotFoundError if no matching file is found.
    """
    sub_str = f"sub-{sub:03d}"
    sub_int = f"{sub:03d}"           # bare zero-padded ID used in filenames

    # Layout: <da_root>/<exp_out>/preproc/<sub_str>/<stage_dir>/
    preproc_root = os.path.join(da_root, exp_out, "preproc", sub_str, stage_dir)

    # Primary: prefix + bare subject integer (e.g. 26BB_62_001_..._base.set)
    pat = os.path.join(preproc_root, f"{out_prefix}*{sub_int}*_base.set")
    hits = glob.glob(pat)

    # Fallback: any file ending in _base.set in the same directory
    if len(hits) == 0 and allow_fallback:
        pat = os.path.join(preproc_root, "*_base.set")
        hits = glob.glob(pat)

    if len(hits) == 0:
        raise FileNotFoundError(
            f"No _base.set found for {sub_str} in {preproc_root}\n"
            f"  Pattern tried: {pat}"
        )

    hits.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return hits[0]

# ====================================================================
# EPOCH LOADING
# ====================================================================
def src_read_epochs(set_path: str) -> mne.Epochs:
    """
    Load an EEGLAB .set file as MNE Epochs with data pre-loaded into memory.

    EEGLAB applies average reference in-place before saving, so the loaded
    data is correctly referenced but MNE sees no reference projector. MNE's
    forward/inverse machinery requires a projector-based reference, so we
    re-register the existing average reference as a projector here. This does
    not alter the data — it just satisfies MNE's bookkeeping requirement.

    Raises RuntimeError if MNE cannot find the reader function.
    """
    # BUG FIX: parameter was named `set_oath` (typo) — corrected to `set_path`
    ep = _read_epochs_eeglab(set_path, verbose="ERROR")
    ep.load_data()
    ep.set_eeg_reference("average", projection=True, verbose="ERROR")
    return ep