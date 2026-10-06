function [feat, pw, f, P] = spec_window_alpha_features(PT, X, timesSec, win, fs, F)
% SPEC_WINDOW_ALPHA_FEATURES  11-metric alpha feature set for one window
% V 1.0.0
%
% Python twin: src_alpha_features.src_compute_window_alpha_features
% (trial x unit part). Per trial x channel:
%   pow_*      mean of the filter-Hilbert power time course inside the window
%   paf_cog_hz CoG of the window's Welch PSD over the alpha band
%   ratios     from the three powers (spec_metrics_from_powers)
%
% Inputs:
%   PT       : struct slow/fast/alpha power time courses [nChan x nTime x nTr]
%   X        : [nChan x nTime x nTr] full-epoch data (cropped here for Welch)
%   timesSec : [1 x nTime] seconds
%   win      : [tmin tmax] seconds
%   fs       : Hz
%   F        : cfg.spectral.feat (bands, psd)
%
% Outputs:
%   feat : struct of [nChan x nTr] metrics (unprefixed)
%   pw   : struct slow/fast/alpha window-mean powers [nChan x nTr]
%   f    : [1 x nF] PSD frequencies
%   P    : [nChan x nF x nTr] per-trial window PSD

m = spec_window_mask(timesSec, win(1), win(2), fs);
if ~any(m)
    error('spec_window_alpha_features:EmptyWindow', ...
        'Window [%.3f %.3f] s has no samples in [%.3f %.3f] s.', win(1), win(2), timesSec(1), timesSec(end));
end

nCh = size(X, 1);
nTr = size(X, 3);
pw = struct();
pw.slow  = reshape(mean(PT.slow(:, m, :), 2),  nCh, nTr);
pw.fast  = reshape(mean(PT.fast(:, m, :), 2),  nCh, nTr);
pw.alpha = reshape(mean(PT.alpha(:, m, :), 2), nCh, nTr);

[f, P] = spec_welch_psd(X(:, m, :), fs, F.psd.fmin_hz, F.psd.fmax_hz, ...
    F.psd.window_sec, F.psd.overlap_frac, F.psd.df_target_hz);
cog = spec_cog_from_psd(f, P, F.bands.alpha_hz);

feat = spec_metrics_from_powers(pw.slow, pw.fast, pw.alpha, cog);
end
