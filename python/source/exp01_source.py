"""
exp01_source.py - Entry point for source localization (Experiment 1)
V 2.0.0

Call Chain:
    exp01_source.py -> source_default.py -> source_core.py
                                         -> src_io / src_assets / src_inverse /
                                            src_spectral / src_prestim /
                                            src_poststim / src_lep /
                                            src_fooof / src_write / src_plot

Usage (command line):
    python exp01_source.py
    python exp01_source.py --subjects 1 2 3 10
    python exp01_source.py --cfg /path/to/exp01.json

Config JSON expected at: <repo_root>/config/exp01.json

Minimum required JSON keys
---------------------------
    cfg.exp.out_prefix                  filename prefix for preprocessed .set files
    cfg.exp.subjects                    list of integer subject IDs
    cfg.source.fsaverage.subjects_dir   path to the directory containing fsaverage/

All other cfg.source keys have safe defaults (see source_default.py).
"""

from __future__ import annotations

import argparse
import json
import os

from source_default import source_default

def exp01_source(subjects_override = None, cfg_path: str = None):
    """
    Entrypoint for Experiment 1 source localization.
    
    Args:
        subjects_override : Optional list of integer subject IDs.
                            When provided, overrides cfg.exp.subjects.
        cfg_path          : Optional explicit path to the experiment JSON config.
                            Defaults to <repo_root>/config/exp01.json.
    """
    exp_id = "exp01"

    if cfg_path is None:
        this_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(os.path.join(this_dir, "..", ".."))
        cfg_path = os.path.join(repo_root, "config", f"{exp_id}.json")

    if not os.path.exists(cfg_path):
        raise FileNotFoundError(
            f"Config JSON not found: {cfg_path}\n"
            "Pass --cfg /path/to/config.json or set cfg_path explicitly."
        )
    
    with open(cfg_path, "r") as fh:
        cfg_in = json.load(fh)

    source_default(exp_id = exp_id, cfg_in = cfg_in, subjects_override = subjects_override)

def _parse_args():
    ap = argparse.ArgumentParser(description = "sLORETA source localization - Experiment 1.")
    ap.add_argument("--cfg", default = None, metavar = "PATH",
                    help = "Path to experiment JSON config (default: <repo>/config/exp01.json)")
    ap.add_argument("--subjects", nargs = "+", type = int, default = None, metavar = "ID",
                    help = "Subject IDs to process (overrides cfg.exp.subjects).")
    return ap.parse_args()

if __name__ == "__main__":
    args = _parse_args()
    exp01_source(subjects_override = args.subjects, cfg_path = args.cfg)