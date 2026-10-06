"""
fooof_bridge.py - MATLAB -> Python FOOOF/specparam bridge (spectral pipeline).
V 2.0.0

Called by matlab/spectral/helpers/spec_run_fooof_python.m. Fits every ROW of
a PSD matrix with python/source/src_fooof.src_fit_fooof, i.e. the exact
function the source pipeline uses, so channel-space and source-space FOOOF
outputs come from one implementation.

V2.0.0 changes vs V1.x:
    - Imports src_fooof from ../source instead of carrying its own copy
      (works with specparam 2.x, specparam 1.x or legacy fooof).
    - Rows are channels (subject x channel trial-mean PSD), not trials.
    - Inputs/outputs are plain numeric CSVs (no JSON NaN handling issues):
        --freq  one row or column of frequencies (Hz)
        --psd   n_rows x n_freqs matrix, no header
        --cfg   JSON {"fooof": {...features.fooof...}, "alpha_band_hz": [8, 12]}
        --out   CSV with header row = FOOOF_NAMES, one line per PSD row

Usage:
    python fooof_bridge.py --freq f.csv --psd psd.csv --cfg cfg.json --out out.csv
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "source")))

from src_fooof import FOOOF_NAMES, fooof_available, fooof_package_name, src_fit_fooof_rows  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", required = True)
    ap.add_argument("--psd", required = True)
    ap.add_argument("--cfg", required = True)
    ap.add_argument("--out", required = True)
    args = ap.parse_args()

    if not fooof_available():
        sys.exit("fooof_bridge: neither 'specparam' nor 'fooof' is installed (pip install specparam).")

    freqs = np.loadtxt(args.freq, delimiter = ",", ndmin = 1).ravel()
    psd = np.loadtxt(args.psd, delimiter = ",", ndmin = 2)
    if psd.shape[1] != freqs.size:
        sys.exit(f"fooof_bridge: psd has {psd.shape[1]} columns but {freqs.size} frequencies.")

    with open(args.cfg, "r") as fh:
        cfg = json.load(fh)
    fooof_cfg = cfg.get("fooof", {})
    alpha_band = tuple(cfg.get("alpha_band_hz", [8.0, 12.0]))

    rows = src_fit_fooof_rows(freqs, psd, fooof_cfg, alpha_band)

    n_fail = 0
    with open(args.out, "w") as fh:
        fh.write(",".join(FOOOF_NAMES) + "\n")
        for r in rows:
            if "fail_reason" in r:
                n_fail += 1
                print(f"[FOOOF] row {r['row']} failed: {r['fail_reason']}", file = sys.stderr)
            fh.write(",".join("NaN" if not np.isfinite(r[k]) else repr(float(r[k])) for k in FOOOF_NAMES) + "\n")

    print(f"[FOOOF] {fooof_package_name()}: fitted {len(rows) - n_fail}/{len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
