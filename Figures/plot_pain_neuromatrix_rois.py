"""
plot_pain_neuromatrix_rois.py - Pain neuromatrix ROI surface render
V 1.0.0

Renders the six cortical pain neuromatrix ROIs on the fsaverage inflated
brain surface using mne.viz.Brain. Produces a 2x2 panel figure (lh lateral,
rh lateral, lh medial, rh medial) with a shared legend mapping colour to ROI.

Output
------
    <OUT_DIR>/pain_neuromatrix_rois.png

Usage
------
    python plot_pain_neuromatrix_rois.py

Edit the CONFIG block below before running.
"""

from __future__ import annotations

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import mne

# =============================================================================
# CONFIG
# =============================================================================

SUBJECTS_DIR = "/home/UWO/darsenea/mne_data/"
OUT_DIR = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/figures"
PARCELLATION = "aparc"

# Pain neuromatrix ROIs - each value is a list of aparc label names
# (both hemispheres included; they will all be coloured the same per ROI)
PAIN_NEUROMATRIX = {
    "S1": ["postcentral-lh", "postcentral-rh",
           "paracentral-lh", "paracentral-rh"],
    "S2": ["supramarginal-lh", "supramarginal-rh"],
    "ACC": ["caudalanteriorcingulate-lh", "caudalanteriorcingulate-rh",
            "rostralanteriorcingulate-lh", "rostralanteriorcingulate-rh"],
    "Insula": ["insula-lh", "insula-rh"],
    "dlPFC": ["rostralmiddlefrontal-lh", "rostralmiddlefrontal-rh",
              "caudalmiddlefrontal-lh", "caudalmiddlefrontal-rh"],
    "M1": ["precentral-lh", "precentral-rh"],
}

# Colourblind-friendly, visually distinct palette (one per ROI)
ROI_COLOURS = {
    "S1": (0.894, 0.102, 0.110), # red
    "S2": (0.216, 0.494, 0.722), # blue
    "ACC": (0.302, 0.686, 0.290), # green
    "Insula": (0.596, 0.306, 0.639), # purple
    "dlPFC": (1.000, 0.498, 0.000), # orange
    "M1": (0.651, 0.337, 0.157), # brown
}

# Brain render settings
SURFACE = "inflated"
ALPHA = 0.9 # ROI opacity
DPI = 300
FIG_SIZE = (16, 10)

# Four views: (hemi, view)
VIEWS = [
    ("lh", "lateral"),
    ("rh", "lateral"),
    ("lh", "medial"),
    ("rh", "medial"),
]

VIEW_TITLES = {
    ("lh", "lateral"): "Left - Lateral",
    ("rh", "lateral"): "Right - Lateral",
    ("lh", "medial"): "Left - Medial",
    ("rh", "medial"): "Right - Medial",
}

# =============================================================================
# HELPERS
# =============================================================================

def _load_target_labels(subjects_dir: str, parc: str, roi_dict: dict) -> dict[str, list]:
    """
    Load mne.Label objects for all target parcels, grouped by ROI name.
    Returns dict mapping roi_name -> list of mne.Label objects.
    """
    all_labels = mne.read_labels_from_annot(
        subject = "fsaverage", parc = parc, subjects_dir = subjects_dir, verbose = "ERROR",
    )
    by_name = {lab.name: lab for lab in all_labels}

    grouped = {}
    for roi_name, label_names in roi_dict.items():
        matched = []
        for name in label_names:
            if name in by_name:
                matched.append(by_name[name])
            else:
                print(f"[WARN] ROI '{roi_name}': label '{name}' not found - skipping")
            if matched:
                grouped[roi_name] = matched
            else:
                print(f"[WARN] ROI '{roi_name}': no labels matched - skipping entirely")

        return grouped
    

def _render_view(subjects_dir: str, hemi: str, view: str, roi_labels: dict[str, list], roi_colors: dict[str, tuple]) -> np.ndarray:
    """
    Render one brain view and return as an RGBA numpy array.
    """
    brain = mne.viz.Brain(
        subject = "fsaverage",
        hemi = hemi,
        surf = SURFACE,
        subjects_dir = subjects_dir,
        background = "black",
        size = (800, 600),
        verbose = "ERROR",
    )

    for roi_name, labels in roi_labels.items():
        color = roi_colors[roi_name]
        for lab in labels:
            if lab.hemi == hemi:
                brain.add_label(lab, color = color, alpha = ALPHA, borders = False)

    brain.show_view(view)
    img = brain.screenshot(time_viewer = False)
    brain.close()
    return img

# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok = True)
    os.environ["SUBJECTS_DIR"] = SUBJECTS_DIR

    print("[LABELS] Loading parcellation labels...")
    roi_labels = _load_target_labels(SUBJECTS_DIR, PARCELLATION, PAIN_NEUROMATRIX)
    print(f"[LABELS] {len(roi_labels)} ROIs loaded: {list(roi_labels.keys())}")

    # -- Render each view ----
    imgs = {}
    for hemi, view in VIEWS:
        key = (hemi, view)
        print(f"[RENDER] {VIEW_TITLES[key]}...")
        imgs[key] = _render_view(SUBJECTS_DIR, hemi, view, roi_labels, ROI_COLOURS)

    # -- Assemble 2x2 panel ----
    fig, axes = plt.subplots(2, 2, figsize = FIG_SIZE, facecolor = "black")
    fig.subplots_adjust(wspace = 0.02, hspace = 0.05, bottom = 0.12)

    panel_order = [
        ("lh", "lateral"), ("rh", "lateral"),
        ("lh", "medial"), ("rh", "medial"),
    ]

    for ax, key in zip(axes.flat, panel_order):
        ax.imshow(imgs[key])
        ax.set_title(VIEW_TITLES[key], color = "white", fontsize = 12, pad = 6)
        ax.axis("off")

    # -- Legend ----
    legend_patches = [
        mpatches.Patch(color = ROI_COLOURS[roi], label = roi)
        for roi in PAIN_NEUROMATRIX.keys()
        if roi in roi_labels
    ]
    fig.legend(
        handles = legend_patches,
        loc = "lower center",
        ncol = len(legend_patches),
        fontsize = 12,
        framealpha = 0.15,
        facecolor = "black",
        edgecolor = "white",
        labelcolor = "white",
        bbox_to_anchor = (0.5, 0.01),
    )

    out_path = os.path.join(OUT_DIR, "pain_neuromatrix_rois.png")
    fig.savefig(out_path, dpi = DPI, bbox_inches = "tight", facecolor = "black")
    plt.close(fig)
    print(f"[FIG] {out_path}")
    print("[DONE]")

if __name__ == "__main__":
    main()