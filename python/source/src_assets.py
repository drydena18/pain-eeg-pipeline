"""
src_assets.py - fsaverage BEM/trans/src asset loading and label/ROI helpers.
V 2.0.0

Covers:
    - Validating and returning the three required fsaverage files
    - Loading all parcellation labels for fsaverage
    - Building macro-ROIs from user-supplied label lists, optionally with
      left- and right-hemisphere versions alongside the bilateral one

V2.0.0 changes vs V1.x:
    - src_build_custom_rois(..., hemisphere_split=True) adds <ROI>_lh and
      <ROI>_rh next to every bilateral macro-ROI, so S1 -> S1, S1_lh, S1_rh.
      The bilateral ROI is unchanged. Each hemisphere ROI is extracted as its
      own label, so with mode = "mean_flip" its sign alignment is computed
      within that hemisphere only; the bilateral time course is therefore
      NOT the average of the two hemisphere time courses.
    - Single-label ROIs are copied before renaming. V1.x renamed the shared
      Label object held in by_name, so e.g. a one-label "Insula_lh" would
      have renamed the parcellation's insula-lh label in place.
    - Duplicate ROI names raise instead of silently producing two columns
      with the same name.
"""

from __future__ import annotations

import os

import mne


# ====================================================================
# FSAVERAGE ASSET LOADING
# ====================================================================
def src_load_fsaverage_assets(subjects_dir: str) -> tuple[str, str, str]:
    """
    Validate and return the fsaverage BEM solution, EEG->MRI transform and
    ico-5 cortical source space under <subjects_dir>/fsaverage/bem/.

    Raises FileNotFoundError listing every missing file.
    """
    bem_dir = os.path.join(subjects_dir, "fsaverage", "bem")
    assets = {
        "BEM solution": os.path.join(bem_dir, "fsaverage-5120-5120-5120-bem-sol.fif"),
        "MRI transform": os.path.join(bem_dir, "fsaverage-trans.fif"),
        "source space": os.path.join(bem_dir, "fsaverage-ico-5-src.fif"),
    }
    missing = [f"[{k}] {v}" for k, v in assets.items() if not os.path.exists(v)]
    if missing:
        raise FileNotFoundError(
            "Missing fsaverage asset(s). Run mne.datasets.fetch_fsaverage() to download them:\n"
            + "\n".join(missing)
        )
    paths = list(assets.values())
    return paths[0], paths[1], paths[2]


# ====================================================================
# LABEL / ROI LOADING
# ====================================================================
def src_load_labels(subjects_dir: str, parcellation: str) -> tuple[list, dict]:
    """
    Load all parcellation labels for fsaverage, dropping 'unknown' parcels.

    Returns:
        labels  : list of mne.Label (both hemispheres)
        by_name : dict label.name -> mne.Label
    """
    labels = mne.read_labels_from_annot(
        subject = "fsaverage",
        parc = parcellation,
        subjects_dir = subjects_dir,
        verbose = "ERROR",
    )
    labels = [lab for lab in labels if not lab.name.startswith("unknown")]
    by_name = {lab.name: lab for lab in labels}
    return labels, by_name


def _norm_label_name(lbl: str) -> str:
    """Normalise 'lh-X' / 'rh-X' to MNE's 'X-lh' / 'X-rh'."""
    s = lbl.strip()
    if s.startswith("lh-"):
        return s[3:] + "-lh"
    if s.startswith("rh-"):
        return s[3:] + "-rh"
    return s


def _combine(labels: list, name: str):
    """Sum labels into one (Label or BiHemiLabel) without mutating inputs."""
    combined = labels[0].copy()
    for lab in labels[1:]:
        combined = combined + lab
    combined.name = name
    return combined


def src_build_custom_rois(
        by_name: dict,
        custom_rois: dict,
        hemisphere_split: bool = False,
) -> tuple[list, list]:
    """
    Merge parcellation labels into macro-ROIs.

    Label names may be written 'postcentral-lh' or 'lh-postcentral'.

    Args:
        by_name          : from src_load_labels()
        custom_rois      : macro-ROI name -> list of label names
        hemisphere_split : if True, every macro-ROI with labels in both
                           hemispheres also yields <name>_lh and <name>_rh,
                           listed immediately after the bilateral ROI

    Returns:
        macro_labels : list of mne.Label / BiHemiLabel
        macro_names  : list of ROI names in the same order

    Unmatched label names are printed as warnings; an ROI with no matched
    labels is skipped.
    """
    macro_labels: list = []
    macro_names: list = []

    def _add(name: str, labs: list):
        if name in macro_names:
            raise ValueError(f"Duplicate ROI name '{name}' (check custom_rois / hemisphere split).")
        macro_labels.append(_combine(labs, name))
        macro_names.append(name)

    for macro_name, parts in custom_rois.items():
        matched, missing = [], []
        for part in parts:
            norm = _norm_label_name(part)
            (matched if norm in by_name else missing).append(norm)
        if missing:
            print(f"[WARN] ROI '{macro_name}': labels not found in parcellation: {missing}")
        if not matched:
            print(f"[WARN] ROI '{macro_name}': no labels matched - skipping entirely.")
            continue

        labs = [by_name[n] for n in matched]
        _add(macro_name, labs)

        if hemisphere_split:
            lh = [lab for lab in labs if lab.hemi == "lh"]
            rh = [lab for lab in labs if lab.hemi == "rh"]
            if lh and rh:
                _add(f"{macro_name}_lh", lh)
                _add(f"{macro_name}_rh", rh)
            else:
                print(f"[INFO] ROI '{macro_name}' is unilateral; no _lh/_rh split added.")

    return macro_labels, macro_names
