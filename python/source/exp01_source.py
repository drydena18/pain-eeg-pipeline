from __future__ import annotations
import os
import argparse
import json

from source_default import source_default

def exp01_source(subjects_override = None, cfg_path = None):
    exp_id = "exp01"

    if cfg_path is None:
        # Adjust to repo layout
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        cfg_path = os.path.join(repo_root, "config", f"{exp_id}.json")

    if not os.path.exists(cfg_path):
        raise FileNotFoundError(f"Config JSON not found: {cfg_path}")
    
    with open(cfg_path, 'r') as f:
        cfg_in = json.load(f)

    # Dispatch
    source_default(exp_id = exp_id, cfg_in = cfg_in, subjects_override = subjects_override)

def _parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfg", default = None, help = "Path to expXX.json")
    ap.add_argument("--subjects", nargs = "+", type = int, default = None)
    return ap.parse_args()

if __name__ == "__main__":
    args = _parse_args()
    exp01_source(subjects_override = args.subjects, cfg_path = args.cfg)