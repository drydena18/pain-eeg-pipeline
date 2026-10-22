"""
src_io.py - IO helpers for the source localization pipeline.
V 2.0.0

Covers:
    - Log file open / write / close
    - Preprocessed .set file resolution
    - EEGLAB .set -> MNE Epochs loading (version-compatible shim)
    - Channel positions for .set files saved without chanlocs, looked up in
      the same BESA .elp preprocessing uses (P.CORE.ELP_FILE)

V2.0.0:
    - src_read_elp_eeglab(): reads the .elp exactly as EEGLAB's
      pop_chanedit(EEG, 'lookup', elp) does (readlocs BESA path +
      convertlocs sphbesa2all / topo2all), giving the same chanlocs X/Y/Z.
    - src_apply_montage(): when the .set has no positions, applies those
      X/Y/Z by label with MNE's own EEGLAB conversion ([-Y, X, Z], mm -> m),
      so positions are identical to every experiment whose 08_base kept its
      chanlocs. Controlled by cfg.source.montage (see source_default.py).
    - src_read_epochs() applies positions BEFORE registering the average-
      reference projector and returns (epochs, montage_tag).

Why this is valid after preprocessing: positions are metadata. Filtering,
notch, resampling, average re-reference, ICA back-projection, epoching and
baseline removal never read them.
"""

from __future__ import annotations

import csv
import datetime
import glob
import os
import re

import numpy as np
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
# CHANNEL POSITIONS (from the same BESA .elp file preprocessing uses)
# ====================================================================
MONTAGE_DEFAULTS = {
    # auto   : keep the .set's positions if present; otherwise look them up in the ELP
    # elp    : always look positions up in the ELP (replaces any in the .set)
    # native : never touch positions; raise if missing
    "mode": "auto",
    "elp_file": "AUTO",             # AUTO -> <da_root>/<exp_out>/resource/standard-10-5-cap385.elp (= P.CORE.ELP_FILE)
    "relabel_csv": "",              # optional raw_label,std_label CSV applied before lookup
    "drop_unmatched": True,         # drop EEG channels the ELP has no entry for
    "max_unmatched_frac": 0.10,     # ...but raise if more than this fraction is unmatched
}

_GENERIC_NAME = re.compile(r"^EEG\s?\d+$", re.IGNORECASE)    # MNE's names when chanlocs are empty
_BIOSEMI_RAW = re.compile(r"^[AB]0*\d{1,2}$", re.IGNORECASE)


def _has_pos(ch) -> bool:
    pos = np.asarray(ch["loc"][:3], dtype = float)
    return bool(np.all(np.isfinite(pos)) and np.linalg.norm(pos) > 1e-6)


def _eeg_channels_without_position(info) -> list[str]:
    """EEG channel names whose location is NaN or all-zero (no chanlocs)."""
    return [info["ch_names"][i] for i in mne.pick_types(info, eeg = True, exclude = [])
            if not _has_pos(info["chs"][i])]


def _read_relabel_csv(path: str) -> dict:
    """raw_label,std_label CSV -> {RAW_LABEL_UPPER: std_label}."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"source.montage.relabel_csv not found: {path}")
    with open(path, newline = "") as fh:
        rows = [{k.strip().lower(): (v or "").strip() for k, v in r.items()} for r in csv.DictReader(fh)]
    if not rows or not {"raw_label", "std_label"} <= set(rows[0]):
        raise ValueError(f"{path} must have columns raw_label,std_label")
    return {r["raw_label"].upper(): r["std_label"] for r in rows if r["raw_label"]}


def src_read_elp_eeglab(elp_path: str) -> dict:
    """
    Read a BESA spherical .elp exactly as EEGLAB's
    pop_chanedit(EEG, 'lookup', elp) does, returning the chanlocs X / Y / Z
    (EEGLAB frame: +X nose, +Y left, +Z up; units of the file radius, mm).

    EEGLAB path replicated (readlocs 'defaultelp','BESA'):
        rows  : type  label  sph_theta_besa  sph_phi_besa  sph_radius
        header: a lone count line ("346") is skipped
        labels: '.' removed
        convertlocs 'sphbesa2all' then 'topo2all':
            sph2topo(method 1) -> topo2sph(method 2) -> sph2topo(method 2)
            -> topo2sph(method 2) -> MATLAB sph2cart(theta, phi, r)

    Returns:
        {label_lower: (label, type, x, y, z)}   (first row wins on duplicates)
    """
    rows = []
    with open(elp_path) as fh:
        for line in fh:
            tok = line.split()
            if len(tok) < 4:                      # header count / blank lines
                continue
            try:
                th_b, ph_b = float(tok[2]), float(tok[3])
            except ValueError:
                continue
            r = float(tok[4]) if len(tok) > 4 else 85.0
            rows.append((tok[1].replace(".", ""), tok[0], th_b, ph_b, r))
    if not rows:
        raise ValueError(f"No BESA spherical rows read from {elp_path}")

    lab = [r[0] for r in rows]
    typ = [r[1] for r in rows]
    th_b = np.array([r[2] for r in rows], float)
    ph_b = np.array([r[3] for r in rows], float)
    rad = np.array([r[4] for r in rows], float)

    # sphbesa2sph: sph2topo([1 th_b ph_b], 1, 1) -> topo2sph([angle radius], 2)
    radius = np.abs(th_b / 180.0) / 1.0
    angle = np.where(th_b >= 0, 90.0 - ph_b, -90.0 - ph_b)
    sph_theta = -angle
    sph_phi = (0.5 - radius) * 180.0
    # sph2all -> sph2topo(method 2): theta / radius (topographic)
    t_theta = -sph_theta
    t_radius = 0.5 - sph_phi / 180.0
    # topo2all -> topo2sph(method 2) recomputes spherical from topographic
    sph_theta = -t_theta
    sph_phi = (0.5 - t_radius) * 180.0
    # sph2cart (MATLAB: azimuth, elevation, r)
    az, el = sph_theta / 180.0 * np.pi, sph_phi / 180.0 * np.pi
    x = rad * np.cos(el) * np.cos(az)
    y = rad * np.cos(el) * np.sin(az)
    z = rad * np.sin(el)

    out = {}
    for i, name in enumerate(lab):
        out.setdefault(name.lower(), (name, typ[i], x[i], y[i], z[i]))
    return out


def _eeglab_xyz_to_mne(xyz: np.ndarray, mean_radius: float) -> np.ndarray:
    """
    mne.io.eeglab: pos = [-Y, X, Z] scaled by montage_units='auto'
    (mean radius >= 25 -> mm, >= 2.5 -> cm, >= 0.25 -> dm, else m).
    """
    scale = 1e-3 if mean_radius >= 25 else 1e-2 if mean_radius >= 2.5 else 1e-1 if mean_radius >= 0.25 else 1.0
    return np.c_[-xyz[:, 1], xyz[:, 0], xyz[:, 2]] * scale


def src_apply_montage(epochs: mne.Epochs, montage_cfg: dict | None = None, logf = None) -> str:
    """
    Make sure every EEG channel has a 3-D position before the forward model.

    Positions are looked up by label in the same BESA .elp that preprocessing
    uses (P.CORE.ELP_FILE), replicating EEGLAB's lookup and MNE's EEGLAB
    reader, so they are identical to what MNE reads for every experiment
    whose 08_base kept its chanlocs. As in EEGLAB's lookup, channels that
    match an ELP row of type EOG become type 'eog' (excluded from the
    forward model, like in the other experiments).

    Returns:
        "native" | "elp"
    """
    mcfg = {**MONTAGE_DEFAULTS, **(montage_cfg or {})}
    mode = str(mcfg["mode"]).lower()
    if mode not in ("auto", "elp", "native"):
        raise ValueError(f"source.montage.mode = '{mcfg['mode']}' (expected auto | elp | native).")

    n_eeg = len(mne.pick_types(epochs.info, eeg = True, exclude = []))
    max_frac = float(mcfg["max_unmatched_frac"])
    missing = _eeg_channels_without_position(epochs.info)

    if mode == "native" or (mode == "auto" and len(missing) <= max_frac * n_eeg):
        if missing and (mode == "native" or not bool(mcfg["drop_unmatched"])):
            raise ValueError(f"{len(missing)}/{n_eeg} EEG channels have no position (e.g. {missing[:5]}).")
        tag = "native"
        src_logmsg(logf, "[MONTAGE] Using channel positions from the .set (%d EEG channels).", n_eeg)
    else:
        if missing:
            src_logmsg(logf, "[MONTAGE] %d/%d EEG channels have no position in the .set (e.g. %s).",
                       len(missing), n_eeg, ", ".join(missing[:5]))

        generic = [ch for ch in epochs.ch_names if _GENERIC_NAME.match(ch)]
        if len(generic) > len(epochs.ch_names) / 2:
            raise ValueError("The .set has no channel labels either (MNE named them 'EEG 000', ...); "
                             "positions cannot be looked up by label.")

        csv_path = str(mcfg.get("relabel_csv") or "").strip()
        if csv_path:
            lut = _read_relabel_csv(csv_path)
            rename = {ch: lut[ch.upper()] for ch in epochs.ch_names if ch.upper() in lut}
            if rename:
                epochs.rename_channels(rename)
            src_logmsg(logf, "[MONTAGE] Relabelled %d channels from %s", len(rename), csv_path)
        raw_biosemi = [ch for ch in epochs.ch_names if _BIOSEMI_RAW.match(ch)]
        if len(raw_biosemi) > 2:
            raise ValueError(f"{len(raw_biosemi)} channels carry raw BioSemi labels "
                             f"({', '.join(raw_biosemi[:6])}, ...); set source.montage.relabel_csv.")

        elp_path = str(mcfg["elp_file"])
        if not os.path.isfile(elp_path):
            raise FileNotFoundError(f"source.montage.elp_file not found: {elp_path}")
        elp = src_read_elp_eeglab(elp_path)

        # Label lookup (case-insensitive, as pop_chanedit 'lookup')
        hits = [(ch, elp[ch.lower()]) for ch in epochs.ch_names if ch.lower() in elp]
        if not hits:
            raise ValueError(f"No channel label of the .set matches {elp_path}.")
        xyz = np.array([[h[2], h[3], h[4]] for _, h in hits], float)
        pos = _eeglab_xyz_to_mne(xyz, float(np.mean(np.linalg.norm(xyz, axis = 1))))

        # ELP type EOG -> MNE 'eog' (pop_chanedit copies the type when empty)
        eog = {ch: "eog" for ch, h in hits if h[1].upper() == "EOG"
               and mne.channel_type(epochs.info, epochs.ch_names.index(ch)) == "eeg"}
        if eog:
            epochs.set_channel_types(eog, verbose = "ERROR")
            src_logmsg(logf, "[MONTAGE] Typed as EOG from the ELP (not used by the forward model): %s",
                       ", ".join(eog))

        epochs.set_montage(mne.channels.make_dig_montage(
            ch_pos = {ch: p for (ch, _), p in zip(hits, pos)}, coord_frame = "head"),
            match_case = False, on_missing = "ignore", verbose = "ERROR")
        n_set = len(mne.pick_types(epochs.info, eeg = True, exclude = [])) - \
            len(_eeg_channels_without_position(epochs.info))
        src_logmsg(logf, "[MONTAGE] ELP positions applied to %d EEG channels from %s.", n_set, elp_path)
        tag = "elp"

    # Whatever is still unpositioned (no ELP entry) is dropped or raises
    n_eeg = len(mne.pick_types(epochs.info, eeg = True, exclude = []))
    still = _eeg_channels_without_position(epochs.info)
    if still:
        if not bool(mcfg["drop_unmatched"]) or len(still) > max_frac * n_eeg:
            raise ValueError(f"{len(still)}/{n_eeg} EEG channels have no position "
                             f"({', '.join(still)}). Check labels / relabel_csv / max_unmatched_frac.")
        epochs.drop_channels(still)
        src_logmsg(logf, "[MONTAGE][WARN] Dropped %d EEG channel(s) with no position: %s",
                   len(still), ", ".join(still))
    return tag


# ====================================================================
# EPOCH LOADING
# ====================================================================
def src_read_epochs(set_path: str, montage_cfg: dict | None = None, logf = None) -> tuple[mne.Epochs, str]:
    """
    Load an EEGLAB .set file as MNE Epochs with data pre-loaded into memory,
    make sure every EEG channel has a position (src_apply_montage), then
    register the average reference as a projector.

    EEGLAB applies average reference in-place before saving, so the loaded
    data is correctly referenced but MNE sees no reference projector. MNE's
    forward/inverse machinery requires a projector-based reference, so we
    re-register the existing average reference as a projector here. This does
    not alter the data — it just satisfies MNE's bookkeeping requirement.
    Positions are applied first so any channel dropped for lacking one never
    enters the projector.

    Returns:
        epochs, montage_tag ("native" | "elp")

    Raises RuntimeError if MNE cannot find the reader function.
    """
    # BUG FIX: parameter was named `set_oath` (typo) — corrected to `set_path`
    ep = _read_epochs_eeglab(set_path, verbose="ERROR")
    ep.load_data()
    montage_tag = src_apply_montage(ep, montage_cfg, logf)
    ep.set_eeg_reference("average", projection=True, verbose="ERROR")
    return ep, montage_tag