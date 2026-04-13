function spec_plot_trial_spectrogram(outPath, GA_S, f, tMs, alpha, subjid, trialIdx)
% SPEC_PLOT_TRIAL_SPECTROGRAM
% Render a channel-averaged spectrogram for one trial.
%
% Inputs:
%   outPath    : full path for the saved PNG
%   GA_S       : [nFreq x nFrames]  linear power (averaged across channels)
%   f          : [nFreq x 1]        frequency vector (Hz)
%   tMs        : [1 x nFrames]      time vector (ms, epoch-referenced; 0 = stimulus)
%   alpha      : cfg.spectral.alpha  struct with .alpha_hz / .slow_hz / .fast_hz
%   subjid     : numeric subject id
%   trialIdx   : numeric trial index

% ---- Convert to dB -------------------------------------------------------
eps0  = 1e-12;
S_dB  = 10 * log10(GA_S + eps0);   % [nFreq x nFrames]

% Robust colour limits: 2nd / 98th percentile across the image so
% outlier frames (edge artefacts at window boundaries) don't crush the
% range.
flat = sort(S_dB(:));
n = numel(flat);
clo = flat(max(1, round(0.02 * n)));
chi = flat(max(n, round(0.98 * n)));
if clo >= chi, chi = clo + 1; end % guard against flat/all-NaN image
clim = [clo chi];

% ---- Figure layout -------------------------------------------------------
h  = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 960 440]);
ax = axes(h, 'Units', 'normalized', 'Position', [0.08 0.12 0.78 0.76]);

% ---- Draw spectrogram ---------------------------------------------------
% imagesc with explicit x/y vectors maps tMs → X axis, f → Y axis.
% 'axis xy' flips Y so low frequency is at the bottom (default puts
% row 1 at top, which would put high Hz at the bottom).
% Colour limits are applied via caxis rather than the 5-arg imasesc form,
% which is not accepted in all MATLAB version when an axes handle is given
imagesc(ax, tMs, f, S_dB);
axis(ax, 'xy');
caxis(ax, clim);
colormap(ax, 'parula');

% Colorbar
cb = colorbar(ax, 'Location', 'eastoutside');
ylabel(cb, 'Power (dB)', 'FontSize', 9);

% Axis labels
xlabel(ax, 'Time (ms)',       'FontSize', 10);
ylabel(ax, 'Frequency (Hz)', 'FontSize', 10);
title(ax, sprintf('sub-%03d  trial-%03d   GA spectrogram (channels averaged)', ...
    subjid, trialIdx), 'Interpreter', 'none', 'FontSize', 10);

% ---- Overlay markers ----------------------------------------------------
% Stimulus onset (t = 0 ms)
if min(tMs) < 0 && max(tMs) > 0
    xline(ax, 0, 'w-',  'LineWidth', 1.5);
end

% Alpha sub-band boundaries (dashed white lines)
% Total alpha band edges
yline(ax, alpha.alpha_hz(1), 'w--', 'LineWidth', 1.0);
yline(ax, alpha.alpha_hz(2), 'w--', 'LineWidth', 1.0);

% Slow/fast boundary (10 Hz) shown as dotted
if isfield(alpha, 'slow_hz') && isfield(alpha, 'fast_hz')
    boundary = alpha.slow_hz(2); % == alpha.fast_hz(1), i.e. 10 Hz
    yline(ax, boundary, 'w:', 'LineWidth', 1.0);
end

% ---- Axes limits ---------------------------------------------------------
xlim(ax, [tMs(1) tMs(end)]);
ylim(ax, [f(1)   f(end)]);

% ---- Export --------------------------------------------------------------
exportgraphics(h, outPath, 'Resolution', 200);
close(h);

end