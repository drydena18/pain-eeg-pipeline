"""
src_render.py - Offscreen 3D brain rendering for the source localization pipeline
V 1.0.0

Produces per-subject 3D surface renders using MNE's Brain object backed by
PyVista in offscreen mode. No display server required.

Two classes of render
----------------------
1. STC time-windowed maps (src_render_stc_timepoints)
    The GA sLORETA source estimate is rendered at three time points:
        - Pre-stimulus mean (average activation over the pre-stim window)
        - N2 peak latency (from the GA LEP, per subject)
        - P2 peak latency (from the GA LEP, per subject)
    Each render is a 4-view layout: left lateral, right lateral, lieft medial
    right medial

2. ROI scalar maps (src_render_roi_scalar)
    A metric scalar (e.g., BI_pre, CoG_pre, delta_ERD) is painted onto the
    cortical surface by filling each ROI's vertices with its GA value.
    Useful for spatially visualizing where a metric is strongest across the
    parcellation. Rendered for a user-specified list of metric columns.

Offscreen rendering
--------------------
PyVista is set to OFF_SCREEN before any Brain object is created. On HPC
(high-performance computing) systems without a GPU, set
cfg.source.render.use_mesa = true in the JSON to force the Mesa software
rendered via the PYOPENGL_PLATFORM = osmesa environment variable

Output filenames
-----------------
    sub-XXX_render_stc_prestim_mean.png
    sub-XXX_render_stc_n2_peak.png
    sub-XXX_render_stc_p2_peak.png
    sub-XXX_render_roi_{metric}.png

Dependencies
-------------
    mne >= 1.0
    pyvista
    pyvistaqt (may be required for some platforms)
    imageio (for screenshot compositing, optional)
"""

from __future__ import annotations

import os
from typing import Optional

import numpy as np

from src_io import src_logmsg

# ====================================================================
# PYVISTA OFFSCREEN SETUP
# ====================================================================
def _configure_offscreen(use_mesa: bool = False):
    """
    Enable PyVista offscreen rendering.

    On headless servers, set use_mesa = True to force the Mesa OSMesa software
    rasterizer via PYOPENGL_PLATFORM=osmesa. This is slower but requires no
    GPU and no X server.
    """
    if use_mesa:
        os.environ["PYOPENGL_PLATFORM"] = "osmesa"
        os.environ["DISPLAY"] = "" # suppress any stray display references

    try:
        import pyvista
        pyvista.OFF_SCREEN = True
        pyvista.start_xvfb() # no-op if xvfb-run is not available
    except Exception:
        pass # pyvista may not be installed; errors surface later at Brain()

# ====================================================================
# VIEW LAYOUT HELPER
# ====================================================================
def _multi_view_screenshot(brain, views: list[tuple[str, str]]) -> np.ndarray:
    """
    Capture multiple views of a Brain object and stitch them into a single
    image row.

    Args:
        brain : mne.viz.Brain instance (already showing the desired data)
        views : List of (hemi, view) tyuples, e.g.
                [("lh", "lat"), ("rh", "lat"), ("lh", "med"), ("rh", "med")]

    Returns:
        img : np.ndarray of shape (H, W*N, 4) RGBA
    """
    shots = []
    for hemi, view in views:
        try:
            brain.show_view(view, hemi = hemi)
            shots.append(brain.screenshot(mode = "rgba", time_viewer = False))
        except Exception:
            pass

    if not shots:
        return np.zeros((400, 400, 4), dtype = np.uint8)
    
    # Pad all shots to the same height before hstack
    max_h = max(s.shape[0] for s in shots)
    padded = []
    for s in shots:
        if s.shape[0] < max_h:
            pad = np.zeros((max_h - s.shape[0], s.shape[1], s.shape[2]), dtype = s.dtype)
            s = np.vstack([s, pad])
        padded.append(s)
    return np.hstack(padded)

def _save_screenshot(img: np.ndarray, path: str):
    """
    Save a numpy RGBA image array to a PNG file.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize = (img.shape[1] / 150, img.shape[0] / 150), dpi = 150)
        ax.imshow(img)
        ax.axis("off")
        fig.tight_layout(pad = 0)
        fig.savefig(path, dip = 150, bbox_inches = "tight")
        plt.close(fig)
    except Exception as e:
        raise RuntimeError(f"Screenshot save failed: {e}") from e
    
# Default 4-view layout used by all renders
_VIEWS_4 = [("lh", "lateral"), ("rh", "lateral"), ("lh", "medial"), ("rh", "medial")]

# ====================================================================
# STC TIME-WINDOWED RENDERS
# ====================================================================
def src_render_stc_timepoints(
        stc_ga,
        ga_lep_rows: list[dict],
        roi_names: list[str],
        times: np.ndarray,
        pre_tmin: float,
        pre_tmax: float,
        subjects_dir: str,
        fig_dir: str,
        sub_str: str,
        clim_pct: tuple[float, float],
        use_mesa: bool,
        logf,
):
    """
    Render the GA sLORETA source estimate at three time points:
        1. Pre-stimulus mean : average of the STC over [pre_tmin, pre_tmax]
        2. N2 peak latency : from ga_lep_rows (first ROI with a valid value)
        3. P2 peak latency : from ga_lep_rows (first ROI with a valid value)

    Each render is a 4-view composite (lh lat / rh lat / lh med / rh med)

    Args:
        stc_ga       : MNE SourceTimeCourse — grand-average evoked STC.
        ga_lep_rows  : Output of src_compute_lep_ga(), for N2/P2 latencies.
        roi_names    : Ordered list of ROI name strings (for labelling only).
        times        : STC time axis in seconds.
        pre_tmin     : Pre-stim window start (s).
        pre_tmax     : Pre-stim window end (s).
        subjects_dir : Path to the directory containing fsaverage/.
        fig_dir      : Output directory for PNG files.
        sub_str      : Subject string (e.g. "sub-001").
        clim_pct     : (lo_pct, hi_pct) percentile bounds for the colour scale.
                       Default (50, 99) suppresses noise while showing peaks.
        use_mesa     : Force Mesa offscreen renderer (for headless HPC servers).
        logf         : Log file handle.
    """
    _configure_offscreen(use_mesa)

    try:
        import mne
    except ImportError:
        src_logmsg(logf, "[WARN] MNE not available - skipping 3D STC renders.")
        return

    # -- Resolve time points -----
    timepoints: list[tuple[str, float]] = [] # (label, t_sec)

    # 1. Pre-stim mean: create a pseudo-STC with mean activation
    pre_mask = (times >= pre_tmin) & (times <= pre_tmax)
    if np.any(pre_mask):
        stc_prestim_mean = stc_ga.copy()
        mean_data = np.mean(stc_ga.data[:, pre_mask], axis = 1, keepdims = True)
        stc_prestim_mean._data = mean_data
        stc_prestim_mean.tmin = 0.0
        stc_prestim_mean.tstep = 1.0
        timepoints.append(("prestim_mean", stc_prestim_mean, None))

    # 2 & 3. N2 and P2 latencies from GA LEP
    def _first_valid(rows: list[dict], key: str) -> Optional[float]:
        """
        Return the first non-NaN value of 'key' across ROI rows (in ms)
        """
        for row in rows:
            v = row.get(key, float("nan"))
            if not np.isnan(float(v)):
                return float(v)
        return None
    
    n2_lat_ms = _first_valid(ga_lep_rows, "n2_lat_ms")
    p2_lat_ms = _first_valid(ga_lep_rows, "p2_lat_ms")

    for label, lat_ms in [("n2_peak", n2_lat_ms), ("p2_peak", p2_lat_ms)]:
        if lat_ms is not None:
            timepoints.append((label, stc_ga, lat_ms / 1000.0))

    if not timepoints:
        src_logmsg(logf, "[WARN] No valid time points for STC renders - skipping.")
        return
    
    # -- Colour scale from the full GA STC -----
    data_vals = np.abs(stc_ga.data)
    clim_lo = float(np.nanpercentile(data_vals, clim_pct[0]))
    clim_hi = float(np.nanpercentile(data_vals, clim_pct[1]))
    clim = {"kind": "value", "lims": [clim_lo, (clim_lo + clim_hi) / 2, clim_hi]}

    # -- Render each time point -----
    for entry in timepoints:
        if len(entry) == 3:
            label, stc_to_plot, t_sec = entry
        else:
            continue

        out_png = os.path.join(fig_dir, f"{sub_str}_render_stc_{label}.png")

        try:
            kwargs = dict(
                subjects = "fsaverage",
                subjects_dir = subjects_dir,
                hemi = "split",
                views = "lateral",
                colormap = "hot",
                clim = clim,
                time_viewer = False,
                show = False,
                backend = "pyvistaqt",
            )
            if t_sec is not None:
                kwargs["initial_time"] = t_sec

            brain = stc_to_plot.plot(**kwargs)
            img = _multi_view_screenshot(brain, _VIEWS_4)
            brain.close()

            _save_screenshot(img, out_png)
            src_logmsg(logf, "[3D] %s", out_png)

        except Exception as e:
            src_logmsg(logf, "[WARN] STC render failed (%s): %s", label, str(e))

# ====================================================================
# ROI SCALAR MAP RENDERS
# ====================================================================
def src_render_roi_scalar(
        ga_rows: list[dict],
        roi_names: list[str],
        labels: list,
        metric_cols: list[str],
        subjects_dir: str,
        fig_dir: str,
        sub_str: str,
        use_mesa: bool,
        logf,
):
    """
    Paint a per-ROI scalar metric onto the fsaverage cortical surface and
    render a 3D 4-view composite PNG for each metric.

    The value for each ROI is looked up from ga_rows, then the vertices
    belonging to that label are filled with the scalar value. Vartices not
    covered by any label (e.g., the unknown/corpus-callosum parcels) are set
    to NaN and rendered transparent

    Args:
        ga_rows     : Grand-average rows from src_write_ga_csv — list of dicts
                      with "roi" and metric columns.
        roi_names   : Ordered list of ROI names (matches labels order).
        labels      : List of mne.Label objects from src_load_labels / custom ROIs.
        metric_cols : Which metric columns to render (e.g.
                      ["BI_pre", "CoG_pre", "delta_ERD", "n2p2_amp"]).
        subjects_dir : Path containing the fsaverage/ folder.
        fig_dir     : Output directory.
        sub_str     : Subject string (e.g. "sub-001").
        use_mesa    : Force Mesa offscreen renderer.
        logf        : Log file handle.
    """
    _configure_offscreen(use_mesa)

    try:
        import mne
    except ImportError:
        src_logmsg(logf, "[WARN] MNE not available - skipping 3D ROI renders.")
        return
    
    # Build roi_name -> metric_value lookup per metric -----
    # ga_rows is a list of dicts; convert to a simple name-keyed dict
    roi_lookup: dict[str, dict] = {row["roi"]: row for row in ga_rows if "roi" in row}

    # -- Load fsaverage source space to get vertex counts -----
    try:
        src = mne.read_source_spaces(
            os.path.join(subjects_dir, "fsaverage", "bem", "fsaverage-ico-5-src.fif"),
            verbose = "ERROR",
        )
    except Exception as e:
        src_logmsg(logf, "[WARN] Cannot load source space for ROI renders: %s", str(e))
        return
    
    n_verts_lh = src[0]["nuse"] # left hemisphere vertex count
    n_verts_rh = src[1]["nuse"] # right hemisphere vertex count
    n_verts = n_verts_lh + n_verts_rh

    for metric in metrics_cols:
        out_png = os.path.join(fig_dir, f"{sub_str}_render_roi_{metric}.png")

        try:
            # -- Build per-vertex data array -----
            vtx_data = np.full(n_verts, np.nan)

            for roi, row in roi_lookup.items():
                val = row.get(metric, float("nan"))
                if np.isnan(float(val)):
                    continue

                # Find the matching label
                matched_label = next(
                    (lab for lab in labels if lab.name == roi), None
                )
                if matched_label in None:
                    continue

                # Map label vertices onto the flat vertex array
                # lh vertices : indices [0, n_verths_lh]
                # rh vertices : indices [n_verts_lh, n_verts_lh + n_verts_rh]
                if matched_label.hemi == "lh":
                    lh_verts = src[0]["vertno"]
                    label_idx = np.intersec1d(matched_label.vertices, lh_verts,
                                              return_indices = True)[2]
                    vtx_data[label_idx] = float(val)
                elif matched_label.hemi == "rh":
                    rh_verts = src[1]["vertno"]
                    label_idx = np.intersect1d(matched_label.vertices, rh_verts,
                                               return_indices = True)[2]
                    vtx_data[n_verts_lh + label_idx] = float(val)

            # Skip if no ROI had a valid value for this metric
            if np.all(np.isnan(vtx_data)):
                src_logmsg(logf, "[WARN] No valid values for metric '%s' - skipping.", metric)
                continue

            # -- Build a fake STC to leverage Brain's colour mapping -----
            # Replace NaN with 0 for rendering; use transparent colormap for 0
            vtx_data_plot = np.nan_to_num(vtx_data, nan = 0.0)

            stc_roi = mne.SourceEstimate(
                data = vtx_data_plot[:, np.newaxis],
                vertices = [src[0]["vertno"], src[1]["vertno"]],
                tmin = 0.0,
                tstep = 1.0,
                subject = "fsaverage",
            )

            # Symmetric colour sclae centred at 0 for signed metrics:
            # one-sided for non-negative metrics (power, amplitude)
            abs_max = float(np.nanmax(np.abs(vtx_data[~np.isnan(vtx_data)])))
            all_positive = float(np.nanmin(vtx_data[~np.isnan(vtx_data)])) >= 0.0

            if all_positive:
                clim = {"kind": "value", "lims": [0, abs_max / 2, abs_max]}
                cmap = "YlOrRd"
            else:
                clim = {"kind": "value", "lims": [-abs_max, 0, abs_max]}
                cmap = "RdBu_r"

            brain = stc_roi.plot(
                subject = "fsaverage",
                subjects_dir = subjects_dir,
                hemi = "split",
                views = "lateral",
                colormap = cmap,
                clim = clim,
                time_viewer = False,
                show = False,
                backend = "pyvistaqt",
                title = f"{sub_str} {metric}",
            )
            img = _multi_view_screenshot(brain, _VIEWS_4)
            brain.close()

            _save_screenshot(img, out_png)
            src_logmsg(logf, "[3D] %s", out_png)

        except Exception as e:
            src_logmsg(logf, "[WARN] ROI scalar render failed (%s): %s", metric, str(e))