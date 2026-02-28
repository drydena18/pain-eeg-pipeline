import argparse
import json
import numpy as np

from fooof import FOOOF

def load_cfg(path):
    with open(path, "r") as f:
        return json.load(f)
    
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", required = True, help = "CSV with freqd (n_freq x 1)")
    ap.add_argument("--psd", required = True, help = "CSV with PSD (n_trials x n_freq), GA per trial")
    ap.add_argument("--cfg", required = True, help = "JSON config for fooof params")
    ap.add_argument("--out", required = True, help = "Output JSON")
    args = ap.parse_args()

    freqs = np.genfromtxt(args.freq, delimiter = ",", skip_header = 1)
    if freqs.ndim > 1:
        freqs = freqs.squeeze()
    
    psd = np.genfromtxt(args.psd, delimiter = ",", skip_header = 1)
    if psd.ndim == 1:
        psd.psd.reshape(1, -1)
    if psd.shape[1] == (len(freqs) + 1):
        psd = psd[:, 1:]

    cfg = load_cfg(args.cfg)

    fooof_cfg = cfg

    fmin = float(fooof_cfg["fmin_hz"])
    fmax = float(fooof_cfg["fmax_hz"])
    pwl = fooof_cfg.get("peak_width_limits", [1.0, 12.0])
    max_peaks = int(fooof_cfg.get("max_n_peaks", 6))
    min_peak_height = float(fooof_cfg.get("min_peak_height", 0.0))
    aperiodic_mode = fooof_cfg.get("aperiodic_mode", "fixed")
    alpha_band = fooof_cfg.get("alpha_band_hz", [8.0, 12.0])
    peak_threshold = float(fooof_cfg.get("peak_threshold", 1.0))
    verbose = bool(fooof_cfg.get("verbose", False))

    out = {
        "summary": {
            "fmin_hz": fmin,
            "fmax_hz": fmax,
            "peak_width_limits": pwl,
            "max_n_peaks": max_peaks,
            "min_peak_height": min_peak_height,
            "peak_threshold": peak_threshold,
            "aperiodic_mode": aperiodic_mode,
            "alpha_band_hz": alpha_band,
            "n_trials": int(psd.shape[0]),
        },
        "trials": []
    }

    # Fit trial-by-trial
    for t in range(psd.shape[0]):
        fm = FOOOF(
            peak_width_limits = pwl,
            max_n_peaks = max_peaks,
            min_peak_height = min_peak_height,
            peak_threshold = peak_threshold,
            aperiodic_mode = aperiodic_mode,
            verbose = verbose
        )

        try:
            fm.fit(freqs, psd[t, :], [fmin, fmax])

            exponent = float(fm.aperiodic_params_[1]) if fm.aperiodic_params_.size >= 2 else np.nan
            offset = float(fm.aperiodic_params_[0]) if fm.aperiodic_params_.size >= 1 else np.nan
            r2 = float(fm.r_squared_)
            err = float(fm.error_)

            # Find strongest alpha peak (if any) among peak params
            # peak_params_: [cf, pw, bw] per peak
            alpha_cf = np.nan
            alpha_pw = np.nan
            alpha_bw = np.nan

            if fm.peak_params_ is not None and len(fm.peak_params_) > 0:
                peaks = np.array(fm.peak_params_)
                # Keep peaks in alpha band
                keep = (peaks[:, 0] >= alpha_band[0]) & (peaks[:, 0] <= alpha_band[1])
                peaks_a = peaks[keep]
                if peaks_a.size > 0:
                    # strongest = max power
                    idx = np.argmax(peaks_a[:, 1])
                    alpha_cf = float(peaks_a[idx, 0])
                    alpha_pw = float(peaks_a[idx, 1])
                    alpha_bw = float(peaks_a[idx, 2])

            out["trials"].append({
                "trial": int(t + 1),
                "aperiodic_offset": offset,
                "aperiodic_exponent": exponent,
                "r2": r2,
                "error": err,
                "alpha_cf": alpha_cf,
                "alpha_pw": alpha_pw,
                "alpha_bw": alpha_bw
            })

        except Exception as e:
            out["trials"].append({
                "trial": int(t + 1),
                "aperiodic_offset": np.nan,
                "aperiodic_exponent": np.nan,
                "r2": np.nan,
                "error": np.nan,
                "alpha_cf": np.nan,
                "alpha_pw": np.nan,
                "alpha_bw": np.nan,
                "fail_reason": str(e)
            })

    with open(args.out, "w") as f:
        json.dump(out, f, indent = 2)

if __name__ == "__main__":
    main()