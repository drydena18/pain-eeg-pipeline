"""
driver_run_all_subjects.py  –  Run run_one_subject.py once per subject, each
in its own fresh subprocess, so the OS fully reclaims memory between subjects.
V 1.0.0

This avoids the OOM "Killed" issue seen when looping over all subjects within
a single long-running Python process. Safe to re-run after a crash — already-
completed subjects are skipped automatically (run_one_subject.py checks for
an existing checkpoint and exits immediately if found).

Usage
-----
    python driver_run_all_subjects.py
"""

from __future__ import annotations

import os
import glob
import subprocess
import sys

DA_ROOT = "/cifs/seminowicz/eegPainDatasets/CNED/da-analysis"

EXPERIMENTS = [
    "26ByBiosemi",
    "142ByBiosemi",
    "29ByANT",
]

WORKER_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_one_subject.py")


def _discover_subjects(da_root: str, exp_name: str) -> list[int]:
    pattern = os.path.join(da_root, exp_name, "source", "sub-*", "fwd", "*_fwd.fif")
    subs = []
    for path in sorted(glob.glob(pattern)):
        sub_str = os.path.basename(path).replace("_fwd.fif", "")
        try:
            subs.append(int(sub_str.replace("sub-", "")))
        except ValueError:
            continue
    return subs


def main():
    total = 0
    ok = 0
    failed = []

    for exp_name in EXPERIMENTS:
        subs = _discover_subjects(DA_ROOT, exp_name)
        print(f"[EXP] {exp_name}: {len(subs)} subjects")

        for sub in subs:
            sub_str = f"sub-{sub:03d}"
            total += 1
            print(f"\n{'='*60}\n[RUN] {exp_name}/{sub_str}\n{'='*60}")

            result = subprocess.run(
                [sys.executable, WORKER_SCRIPT, exp_name, str(sub)],
                capture_output=False,   # let worker's prints stream directly to console
            )

            if result.returncode == 0:
                ok += 1
            else:
                failed.append((exp_name, sub_str, result.returncode))
                print(f"[DRIVER] {exp_name}/{sub_str} exited with code {result.returncode}")

    print(f"\n{'='*60}")
    print(f"[SUMMARY] {ok}/{total} subjects completed successfully")
    if failed:
        print(f"[SUMMARY] {len(failed)} failed:")
        for exp_name, sub_str, code in failed:
            print(f"    {exp_name}/{sub_str}  (exit code {code})")
        print("\nRe-run this driver script to retry failed/incomplete subjects "
              "— completed checkpoints are skipped automatically.")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()