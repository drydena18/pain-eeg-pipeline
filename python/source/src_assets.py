"""
src_assets.py - fsaverage BEM/trans/src asset loading and label/ROI helpers.

Covers:
    - Validating and returning the three required fsaverage files
    - Loading all parcellation labels for fsaverage
    - Building macro-ROIs from user-supplied label lists
"""

from __future__ import annotations

import os

import mne

# ====================================================================
# FSAVERAGE ASSET LOADING
# ====================================================================
def src_load_fsaverage_assets(subjects_dir: str) -> tuple[str, str, str]:
    """
    Validate and return paths to the three fsaverage BEM (boundary element model)
    assets required by MNE for source localization without subject MRI.

    Expected layout under <subjects_dir>/fsaverage/bem/:
        fsaverage-5120-5120-5120-bem-sol.fif (pre-computed BEM solution)
        fsaverage-trans.fif (EEG -> MRI coordinate transform)
        fsaverage-ico-5-src.fif (cortical source space, ico = 5)

    These files are distributed by MNE via mne.datasets.fetch_fsaverage().

    Returns:
        (bem_sol_path, trans_path, src_space_path) - all validated to exist.

    Raises:
        FileNotFoundError listing every missing file.
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
    return paths[0], paths[1], paths[2] # bem_sol, trans, src_space

# ====================================================================
# LABEL / ROI LOADING
# ====================================================================
def src_load_labels(
        subjects_dir: str,
        parcellation: str,
) -> tuple[list, dict]:
    """
    Load all parcellation labels for the fsaverage template brain.

    Args:
        subjects_dir : Directory containing the fsaverage/ folder
        parcellation : Parcellation name (e.g., "aparc" for Desikan-Killiany,
                       "aparc.a2009s" for Destrieux)

    Returns:
        labels : List of mne.Label objects (all regions, both hemispheres)
        by_name : Dict mapping label.name -> mne.Label for fast lookup
    """
    labels = mne.read_labels_from_annot(
        subject = "fsaverage",
        parc = parcellation,
        subjects_dir = subjects_dir,
        verbose = "ERROR",
    )
    by_name = {lab.name: lab for lab in labels}
    return labels, by_name

def src_build_custom_rois(
        by_name: dict,
        custom_rois: dict,
) -> tuple[list, list]:
    """
    Merge individual parcellation labels into macro-ROIs (regions of interest)
    as defined by the user in the JSON config.

    Label names are accepted in both MNE convention ("postcentral-lh") and the
    reversed "hemisphere-first" form ("lh-postcentral") - both are normalized
    automatically to MNE's convention before lookup.

    Example custom_rois JSON:
        {
            "S1": ["postcentral-lh", "postcentral-rh"],
            "M1": ["precentral-lh", "precentral-rh"],
            "ACC": ["caudal-accumbens-area-lh", "caudal-accumbens-area-rh"]
        }

    Args:
        by_name : Dict from src_load_labels() mapping name -> mne.Label
        custom_rois : Dict mapping macro-ROI name -> list of label name strings

    Returns:
        macro_labels : List of combined mne.Label objects (one per macro-ROI)
        macro_names : List of string names in matching order.

    Warnings are printed (not raised) for unrecognized label names so that a
    partially matched ROI still contributes rather than being silently dropped.
    """
    def _norm(lbl: str) -> str:
        """
        Normalize lh-X / rh-X -> X-lh / X-rh.
        """
        s = lbl.strip()
        if s.startswith("lh-"):
            return s[3:] + "-lh"
        if s.startswith("rh-"):
            return s[3:] + "-rh"
        return s
    
    macro_labels: list = []
    macro_names: list = []
    
    for macro_name, parts in custom_rois.items():
        matched, missing = [], []
        for part in parts:
            norm = _norm(part)
            if norm in by_name:
                matched.append(by_name[norm])
            else:
                missing.append(norm)

        if missing:
            print(f"[WARN] ROI '{macro_name}': labels not found in parcellation: {missing}")
        if not matched:
            print(f"[WARN] ROI '{macro_name}': no labels matched - skipping entirely.")
            continue

        combined = matched[0]
        for lab in matched[1:]:
            combined = combined + lab
        combined.name = macro_name

        macro_labels.append(combined)
        macro_names.append(macro_name)

    return macro_labels, macro_names