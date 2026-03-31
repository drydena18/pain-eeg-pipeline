import argparse
import json
import numpy as np

from fooof import FOOOF

def load_cfg(path):
    with open(path, "r") as f:
        return json.load(f)
    
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", required = True, help = "CSV with freqs (n_freq x 1)")
    ap.add_argument("--psd",  required = True, help = "CSV with PSD (n_trials x n_freq), GA per trial")
    ap.add_argument("--cfg",  required = True, help = "JSON config for FOOOF params")
    ap.add_argument("--out",  required = True, help = "Output JSON")
    args = ap.parse_args()

    freqs = np.genfromtxt(args.freq, delimiter = ",", skip_header = 1)
    if freqs.ndim > 1:
        freqs = freqs.squeeze()
    
    psd = np.genfromtxt(args.psd, delimiter = ",", skip_header = 1)
    if psd.ndim == 1:
        psd = psd.reshape(1, -1)
    if psd.shape[1] == (len(freqs) + 1):
        # Drop leading trial-index column if present
        psd = psd[:, 1:]

    cfg = load_cfg(args.cfg)

    fmin             = float(cfg["fmin_hz"])
    fmax             = float(cfg["fmax_hz"])
    pwl              = cfg.get("peak_width_limits", [1.0, 12.0])
    max_peaks        = int(cfg.get("max_n_peaks", 10))
    min_peak_height  = float(cfg.get("min_peak_height", 0.0))
    aperiodic_mode   = cfg.get("aperiodic_mode", "fixed")
    alpha_band       = cfg.get("alpha_band_hz", [8.0, 12.0])
    peak_threshold   = float(cfg.get("peak_threshold", 1.0))
    verbose          = bool(cfg.get("verbose", False))

    out = {
        "summary": {
            "fmin_hz":           fmin,
            "fmax_hz":           fmax,
            "peak_width_limits": pwl,
            "max_n_peaks":       max_peaks,
            "min_peak_height":   min_peak_height,
            "peak_threshold":    peak_threshold,
            "aperiodic_mode":    aperiodic_mode,
            "alpha_band_hz":     alpha_band,
            "n_trials":          int(psd.shape[0]),
        },
        "trials": []
    }

    for t in range(psd.shape[0]):
        fm = FOOOF(
            peak_width_limits = pwl,
            max_n_peaks       = max_peaks,
            min_peak_height   = min_peak_height,
            peak_threshold    = peak_threshold,
            aperiodic_mode    = aperiodic_mode,
            verbose           = verbose
        )

        try:
            fm.fit(freqs, psd[t, :], [fmin, fmax])

            # ---- Aperiodic parameters ----
            # fixed mode: [offset, exponent]
            # knee  mode: [offset, knee, exponent]
            ap_params = fm.aperiodic_params_

            offset   = float(ap_params[0]) if ap_params.size >= 1 else float("nan")
            knee_val = float("nan")
            exponent = float("nan")

            if aperiodic_mode == "knee" and ap_params.size >= 3:
                knee_val = float(ap_params[1])
                exponent = float(ap_params[2])
            elif ap_params.size >= 2:
                exponent = float(ap_params[1])

            r2  = float(fm.r_squared_)
            err = float(fm.error_)

            # ---- Alpha peak (strongest in band) ----
            alpha_cf = float("nan")
            alpha_pw = float("nan")
            alpha_bw = float("nan")

            if fm.peak_params_ is not None and len(fm.peak_params_) > 0:
                peaks = np.array(fm.peak_params_)
                keep  = (peaks[:, 0] >= alpha_band[0]) & (peaks[:, 0] <= alpha_band[1])
                peaks_a = peaks[keep]
                if peaks_a.size > 0:
                    idx      = np.argmax(peaks_a[:, 1])
                    alpha_cf = float(peaks_a[idx, 0])
                    alpha_pw = float(peaks_a[idx, 1])
                    alpha_bw = float(peaks_a[idx, 2])

            out["trials"].append({
                "trial":              int(t + 1),
                "aperiodic_offset":   offset,
                "aperiodic_knee":     knee_val,   # NaN for fixed mode
                "aperiodic_exponent": exponent,
                "r2":                 r2,
                "error":              err,
                "alpha_cf":           alpha_cf,
                "alpha_pw":           alpha_pw,
                "alpha_bw":           alpha_bw,
            })

        except Exception as e:
            out["trials"].append({
                "trial":              int(t + 1),
                "aperiodic_offset":   float("nan"),
                "aperiodic_knee":     float("nan"),
                "aperiodic_exponent": float("nan"),
                "r2":                 float("nan"),
                "error":              float("nan"),
                "alpha_cf":           float("nan"),
                "alpha_pw":           float("nan"),
                "alpha_bw":           float("nan"),
                "fail_reason":        str(e),
            })

    with open(args.out, "w") as f:
        json.dump(out, f, indent = 2)

if __name__ == "__main__":
    main()