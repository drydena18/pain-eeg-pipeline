"""
plot_pain_neuromatrix_rois.py  -  Pain neuromatrix ROI surface render.
V 2.0.0

Uses PyVista directly for offscreen rendering (no Qt/display required).
Loads the fsaverage inflated surface mesh and colors ROI vertices from
the aparc annotation file.

Output
------
    <OUT_DIR>/pain_neuromatrix_rois.png

Usage
-----
    python plot_pain_neuromatrix_rois.py
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import nibabel as nib
import pyvista as pv

pv.global_theme.multi_rendering_splitting_position = None
os.environ["DISPLAY"] = os.environ.get("DISPLAY", ":99")   # headless fallback

# =============================================================================
# CONFIG
# =============================================================================

SUBJECTS_DIR = "/home/UWO/darsenea/mne_data/"
FSAVERAGE    = os.path.join(SUBJECTS_DIR, "fsaverage")
OUT_DIR      = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis/figures"

PAIN_NEUROMATRIX = {
    "S1":     ["postcentral-lh", "postcentral-rh",
               "paracentral-lh", "paracentral-rh"],
    "S2":     ["supramarginal-lh", "supramarginal-rh"],
    "ACC":    ["caudalanteriorcingulate-lh", "caudalanteriorcingulate-rh",
               "rostralanteriorcingulate-lh", "rostralanteriorcingulate-rh"],
    "Insula": ["insula-lh", "insula-rh"],
    "dlPFC":  ["rostralmiddlefrontal-lh", "rostralmiddlefrontal-rh",
               "caudalmiddlefrontal-lh", "caudalmiddlefrontal-rh"],
    "M1":     ["precentral-lh", "precentral-rh"],
}

ROI_COLORS = {
    "S1":     (0.894, 0.102, 0.110),
    "S2":     (0.216, 0.494, 0.722),
    "ACC":    (0.302, 0.686, 0.290),
    "Insula": (0.596, 0.306, 0.639),
    "dlPFC":  (1.000, 0.498, 0.000),
    "M1":     (0.651, 0.337, 0.157),
}

BRAIN_COLOR  = (0.78, 0.78, 0.78)   # light grey for unlabeled cortex
BACKGROUND   = (0.1, 0.1, 0.1)      # dark background
WINDOW_SIZE  = (1200, 900)
DPI          = 300
FIG_SIZE     = (16, 10)

# Camera positions for each view (hemi, view_name, camera_position)
# PyVista camera: (position, focal_point, view_up)
VIEWS = {
    ("lh", "lateral"): dict(
        position=(-500, 0, 0), focal=(0, 0, 0), up=(0, 0, 1)
    ),
    ("rh", "lateral"): dict(
        position=(500, 0, 0), focal=(0, 0, 0), up=(0, 0, 1)
    ),
    ("lh", "medial"): dict(
        position=(500, 0, 0), focal=(0, 0, 0), up=(0, 0, 1)
    ),
    ("rh", "medial"): dict(
        position=(-500, 0, 0), focal=(0, 0, 0), up=(0, 0, 1)
    ),
}

VIEW_TITLES = {
    ("lh", "lateral"): "Left — Lateral",
    ("rh", "lateral"): "Right — Lateral",
    ("lh", "medial"):  "Left — Medial",
    ("rh", "medial"):  "Right — Medial",
}

PANEL_ORDER = [
    ("lh", "lateral"), ("rh", "lateral"),
    ("lh", "medial"),  ("rh", "medial"),
]

# =============================================================================
# HELPERS
# =============================================================================

def _load_surface(fsaverage: str, hemi: str) -> tuple[np.ndarray, np.ndarray]:
    """Load inflated surface coords and faces from FreeSurfer binary."""
    surf_path = os.path.join(fsaverage, "surf", f"{hemi}.inflated")
    coords, faces = nib.freesurfer.read_geometry(surf_path)
    return coords, faces


def _load_annot(fsaverage: str, hemi: str, parc: str = "aparc"):
    """Load parcellation annotation — returns (labels, ctab, names)."""
    annot_path = os.path.join(fsaverage, "label", f"{hemi}.{parc}.annot")
    labels, ctab, names = nib.freesurfer.read_annot(annot_path)
    names = [n.decode() if isinstance(n, bytes) else n for n in names]
    return labels, ctab, names


def _build_vertex_colors(vertex_labels: np.ndarray, names: list[str],
                          hemi: str, roi_dict: dict, roi_colors: dict,
                          brain_color: tuple) -> np.ndarray:
    """
    Build per-vertex RGB color array.
    ROI vertices get their assigned color; everything else gets brain_color.
    """
    n_verts = len(vertex_labels)
    colors = np.tile(np.array(brain_color, dtype=np.float32), (n_verts, 1))

    for roi_name, label_names in roi_dict.items():
        color = np.array(roi_colors[roi_name], dtype=np.float32)
        for label_name in label_names:
            # Strip hemisphere suffix to get bare parcel name
            bare = label_name.replace(f"-{hemi}", "")
            if f"{bare}-{hemi}" not in label_name and bare not in label_name:
                continue
            # Find index of this parcel in the annotation names
            parcel_name = bare  # e.g. "postcentral"
            if parcel_name in names:
                idx = names.index(parcel_name)
                mask = vertex_labels == idx
                colors[mask] = color

    return colors


def _render_hemi_view(fsaverage: str, hemi: str, view_key: tuple,
                       roi_dict: dict, roi_colors: dict) -> np.ndarray:
    """Render one hemisphere from one camera angle, return screenshot array."""
    coords, faces = _load_surface(fsaverage, hemi)
    vertex_labels, ctab, names = _load_annot(fsaverage, hemi)

    vertex_colors = _build_vertex_colors(
        vertex_labels, names, hemi, roi_dict, roi_colors, BRAIN_COLOR
    )

    # Build PyVista mesh
    n_faces = faces.shape[0]
    pv_faces = np.hstack([np.full((n_faces, 1), 3, dtype=np.int_), faces])
    mesh = pv.PolyData(coords, pv_faces)
    mesh["colors"] = (vertex_colors * 255).astype(np.uint8)

    plotter = pv.Plotter(off_screen=True, window_size=WINDOW_SIZE)
    plotter.background_color = BACKGROUND
    plotter.add_mesh(mesh, scalars="colors", rgb=True, smooth_shading=True)

    cam = VIEWS[view_key]
    plotter.camera.position = cam["position"]
    plotter.camera.focal_point = cam["focal"]
    plotter.camera.up = cam["up"]
    plotter.reset_camera()

    img = plotter.screenshot(return_img=True)
    plotter.close()
    return img


# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("[RENDER] Rendering 4 views...")
    imgs = {}
    for view_key in PANEL_ORDER:
        hemi, view_name = view_key
        print(f"  {VIEW_TITLES[view_key]}...")
        imgs[view_key] = _render_hemi_view(
            FSAVERAGE, hemi, view_key, PAIN_NEUROMATRIX, ROI_COLORS
        )

    # ── Assemble 2x2 panel ────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=FIG_SIZE,
                              facecolor=BACKGROUND)
    fig.subplots_adjust(wspace=0.02, hspace=0.05, bottom=0.12)

    for ax, key in zip(axes.flat, PANEL_ORDER):
        ax.imshow(imgs[key])
        ax.set_title(VIEW_TITLES[key], color="white", fontsize=12, pad=6)
        ax.axis("off")

    # ── Legend ────────────────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color=ROI_COLORS[roi], label=roi)
        for roi in PAIN_NEUROMATRIX
    ]
    fig.legend(
        handles=legend_patches,
        loc="lower center",
        ncol=len(legend_patches),
        fontsize=12,
        framealpha=0.15,
        facecolor=BACKGROUND,
        edgecolor="white",
        labelcolor="white",
        bbox_to_anchor=(0.5, 0.01),
    )

    out_path = os.path.join(OUT_DIR, "pain_neuromatrix_rois.png")
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight",
                facecolor=BACKGROUND)
    plt.close(fig)
    print(f"[FIG] {out_path}")
    print("[DONE]")

if __name__ == "__main__":
    main()