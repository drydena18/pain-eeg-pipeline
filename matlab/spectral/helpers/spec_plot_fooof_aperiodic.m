function spec_plot_fooof_aperiodic(outPath, fooofOut, subjid, cfg, f, gaPxx)
% SPEC_PLOT_FOOOF_APERIODIC  Dedicated aperiodic-fit figure for FOOOF output
% V 1.0.0
%
% Produces a focused 2-row figure covering everything related to the
% aperiodic component of the FOOOF (Fitting Oscillations & One-Over-F) model:
%
%   Row 1 — Scalar time-series (one point per trial):
%     [1,1] Aperiodic Exponent over trials
%     [1,2] Aperiodic Offset over trials
%     [1,3] Aperiodic Knee over trials (knee mode only)
%            — shows R² over trials for fixed mode instead.
%
%   Row 2 — Spectral panels (computed from the mean PSD across trials):
%     [2,1] Mean raw PSD (dB) with the reconstructed trial-averaged
%           aperiodic fit overlaid in red.  Only the FOOOF fitting range
%           [fmin, fmax] is shown.
%     [2,2] Flattened spectrum: raw PSD minus the aperiodic fit (dB).
%           Positive residuals are periodic components (peaks); values
%           near zero confirm aperiodic-dominated regions.  The alpha
%           band [8–12 Hz] is shaded in light blue for orientation.
%     [2,3] Per-trial aperiodic fit residual (root mean squared error —
%           RMSE — of the fit within the fitting range) plotted as a
%           time-series.  Drift upward over the session can indicate
%           progressive changes in the 1/f slope (e.g. fatigue).
%
% Aperiodic reconstruction (from fitted parameters):
%   fixed mode: log10(power_ap) = offset − exponent × log10(f)
%   knee  mode: log10(power_ap) = offset − log10(knee + f^exponent)
%
% Inputs:
%   outPath  : full path for the saved PNG
%   fooofOut : struct from spec_run_fooof_python / spec_fill_fooof_alpha
%   subjid   : integer subject ID
%   cfg      : pipeline cfg struct (optional; used for plot-styling overrides
%              under cfg.spectral.qc.fooof_plot)
%   f        : [1 x nFreq] frequency vector (Hz) — full-range PSD freqs
%              Pass [] to skip row-2 spectral panels.
%   gaPxx    : [nFreq x nTrials] grand-average PSD (linear power, μV²/Hz)
%              Pass [] to skip row-2 spectral panels.

if nargin < 4, cfg   = struct(); end
if nargin < 5, f     = [];       end
if nargin < 6, gaPxx = [];       end

% ---------------------------------------------------------------
% Styling defaults
% ---------------------------------------------------------------
dpi     = 250;
fsTitle = 11;
fsAxis  = 10;
fsTick  = 9;
lw      = 1.25;
mk      = 'o';
ms      = 3.5;

try
    if isfield(cfg, 'spectral') && isfield(cfg.spectral, 'qc') && ...
            isfield(cfg.spectral.qc, 'fooof_plot')
        fp = cfg.spectral.qc.fooof_plot;
        if isfield(fp, 'dpi'),        dpi     = fp.dpi;        end
        if isfield(fp, 'fsTitle'),    fsTitle = fp.fsTitle;    end
        if isfield(fp, 'fsAxis'),     fsAxis  = fp.fsAxis;     end
        if isfield(fp, 'fsTick'),     fsTick  = fp.fsTick;     end
        if isfield(fp, 'lineWidth'),  lw      = fp.lineWidth;  end
        if isfield(fp, 'markerSize'), ms      = fp.markerSize; end
    end
catch
end

% ---------------------------------------------------------------
% Validate input
% ---------------------------------------------------------------
hasFooof = isstruct(fooofOut) && isfield(fooofOut, 'trials') && ~isempty(fooofOut.trials);

h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 1600 800]);
set(h, 'Color', 'w');

if ~hasFooof
    axis off;
    text(0.5, 0.5, sprintf('sub-%03d  FOOOF aperiodic: not available / failed', subjid), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', fsTitle, 'Interpreter', 'none');
    exportgraphics(h, outPath, 'Resolution', dpi);
    close(h);
    return;
end

tr = fooofOut.trials;
nT = numel(tr);

% ---------------------------------------------------------------
% Extract per-trial scalar arrays
% ---------------------------------------------------------------
expnt = nan(nT, 1);
offs  = nan(nT, 1);
knee  = nan(nT, 1);
r2    = nan(nT, 1);
err   = nan(nT, 1);

for t = 1:nT
    if isfield(tr(t), 'aperiodic_exponent'), expnt(t) = tr(t).aperiodic_exponent; end
    if isfield(tr(t), 'aperiodic_offset'),   offs(t)  = tr(t).aperiodic_offset;   end
    if isfield(tr(t), 'aperiodic_knee'),     knee(t)  = tr(t).aperiodic_knee;     end
    if isfield(tr(t), 'r2'),                 r2(t)    = tr(t).r2;                 end
    if isfield(tr(t), 'error'),              err(t)   = tr(t).error;              end
end

tt = (1:nT)';

% ---------------------------------------------------------------
% Determine aperiodic mode and fitting range
% ---------------------------------------------------------------
apMode  = 'fixed';
fmin_fit = NaN;
fmax_fit = NaN;
alphaBand = [8 12];

if isfield(fooofOut, 'summary')
    s = fooofOut.summary;
    if isfield(s, 'aperiodic_mode'), apMode   = char(string(s.aperiodic_mode)); end
    if isfield(s, 'fmin_hz'),        fmin_fit  = s.fmin_hz;                     end
    if isfield(s, 'fmax_hz'),        fmax_fit  = s.fmax_hz;                     end
    if isfield(s, 'alpha_band_hz'),  alphaBand = s.alpha_band_hz;               end
end

useKnee   = strcmpi(apMode, 'knee') && any(~isnan(knee));
mean_exp  = mean(expnt, 'omitnan');
mean_offs = mean(offs,  'omitnan');
mean_knee = mean(knee,  'omitnan');

% ---------------------------------------------------------------
% Layout: 2 rows x 3 cols
% ---------------------------------------------------------------
tl = tiledlayout(h, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

modeLabel = apMode;
if ~isnan(fmin_fit) && ~isnan(fmax_fit)
    modeLabel = sprintf('%s | fit: %.1f–%.1f Hz', apMode, fmin_fit, fmax_fit);
end
title(tl, sprintf('sub-%03d  FOOOF Aperiodic Component (%s)', subjid, modeLabel), ...
    'Interpreter', 'none', 'FontSize', fsTitle + 1, 'FontWeight', 'bold');

% ================================================================
% ROW 1 — Scalar time-series
% ================================================================

% [1,1] Exponent
ax = nexttile(tl, 1);
plot(ax, tt, expnt, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms, 'Color', [0.13 0.47 0.71]);
xlabel(ax, 'Trial',    'FontSize', fsAxis);
ylabel(ax, 'Exponent', 'FontSize', fsAxis);
title(ax, 'Aperiodic Exponent', 'Interpreter', 'none', 'FontSize', fsAxis);
set(ax, 'FontSize', fsTick); grid(ax, 'on');

% Add a smoothed trend line to make session drift visible
if nT >= 5
    hold(ax, 'on');
    sm = smoothdata(expnt, 'movmean', max(3, round(nT / 10)));
    plot(ax, tt, sm, '-', 'LineWidth', lw + 0.5, 'Color', [0.84 0.15 0.16]);
    legend(ax, {'Trial', 'Smoothed'}, 'Location', 'best', 'FontSize', fsTick - 1);
    hold(ax, 'off');
end

% [1,2] Offset
ax = nexttile(tl, 2);
plot(ax, tt, offs, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms, 'Color', [0.17 0.63 0.17]);
xlabel(ax, 'Trial',  'FontSize', fsAxis);
ylabel(ax, 'Offset', 'FontSize', fsAxis);
title(ax, 'Aperiodic Offset', 'Interpreter', 'none', 'FontSize', fsAxis);
set(ax, 'FontSize', fsTick); grid(ax, 'on');

if nT >= 5
    hold(ax, 'on');
    sm = smoothdata(offs, 'movmean', max(3, round(nT / 10)));
    plot(ax, tt, sm, '-', 'LineWidth', lw + 0.5, 'Color', [0.84 0.15 0.16]);
    hold(ax, 'off');
end

% [1,3] Knee (knee mode) or R² (fixed mode)
ax = nexttile(tl, 3);
if useKnee
    plot(ax, tt, knee, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms, 'Color', [0.58 0.40 0.74]);
    xlabel(ax, 'Trial', 'FontSize', fsAxis);
    ylabel(ax, 'Knee',  'FontSize', fsAxis);
    title(ax, 'Aperiodic Knee', 'Interpreter', 'none', 'FontSize', fsAxis);
else
    plot(ax, tt, r2, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms, 'Color', [0.58 0.40 0.74]);
    yline(ax, 0.9, '--k', 'R²=0.9', 'FontSize', fsTick - 1, 'LabelHorizontalAlignment', 'left');
    xlabel(ax, 'Trial', 'FontSize', fsAxis);
    ylabel(ax, 'R²',    'FontSize', fsAxis);
    title(ax, 'Fit R² over Trials', 'Interpreter', 'none', 'FontSize', fsAxis);
    ylim(ax, [max(0, min(r2, [], 'omitnan') - 0.05), 1.02]);
end
set(ax, 'FontSize', fsTick); grid(ax, 'on');

% ================================================================
% ROW 2 — Spectral panels
% ================================================================
hasPSD    = ~isempty(f) && ~isempty(gaPxx) && size(gaPxx, 1) == numel(f);
hasRange  = ~isnan(fmin_fit) && ~isnan(fmax_fit);

if hasPSD && hasRange
    % Mean raw PSD across trials (linear power)
    meanPsd = mean(gaPxx, 2, 'omitnan');   % [nFreq x 1]

    % Restrict to FOOOF fitting range
    idxFit  = (f(:)' >= fmin_fit) & (f(:)' <= fmax_fit);
    f_fit   = f(idxFit);
    psd_fit = meanPsd(idxFit);

    % Reconstruct trial-averaged aperiodic fit
    if useKnee
        log10_ap = mean_offs - log10(mean_knee + f_fit(:) .^ mean_exp);
    else
        log10_ap = mean_offs - mean_exp .* log10(f_fit(:));
    end

    psd_fit_dB = 10 * log10(max(psd_fit(:), 1e-30));
    ap_fit_dB  = 10 * log10(max(10 .^ log10_ap, 1e-30));
    flat_dB    = psd_fit_dB - ap_fit_dB;

    % Per-trial fit RMSE within the fitting range
    trialRMSE = nan(nT, 1);
    for t = 1:nT
        if ~isnan(offs(t)) && ~isnan(expnt(t))
            psd_t = gaPxx(idxFit, t);
            if useKnee && ~isnan(knee(t))
                la_t = offs(t) - log10(knee(t) + f_fit(:) .^ expnt(t));
            else
                la_t = offs(t) - expnt(t) .* log10(f_fit(:));
            end
            diff_dB = 10 * log10(max(psd_t(:), 1e-30)) - 10 * log10(max(10 .^ la_t, 1e-30));
            trialRMSE(t) = sqrt(mean(diff_dB .^ 2, 'omitnan'));
        end
    end

    % [2,1] Raw PSD + aperiodic fit overlay
    ax = nexttile(tl, 4);
    plot(ax, f_fit, psd_fit_dB, 'k-',  'LineWidth', lw,      'DisplayName', 'Raw PSD (mean)');
    hold(ax, 'on');
    plot(ax, f_fit, ap_fit_dB,  'r--', 'LineWidth', lw + 0.5, 'DisplayName', sprintf('Aperiodic fit (%s)', apMode));
    hold(ax, 'off');
    xlabel(ax, 'Frequency (Hz)',         'FontSize', fsAxis);
    ylabel(ax, 'Power (dB re μV²/Hz)',   'FontSize', fsAxis);
    title(ax, 'Mean PSD + Aperiodic Fit', 'Interpreter', 'none', 'FontSize', fsAxis);
    legend(ax, 'Location', 'northeast', 'FontSize', fsTick - 1);
    xlim(ax, [fmin_fit fmax_fit]);
    set(ax, 'FontSize', fsTick); grid(ax, 'on');

    % [2,2] Flattened (periodic) spectrum
    ax = nexttile(tl, 5);
    alphaInRange = alphaBand(1) >= fmin_fit && alphaBand(2) <= fmax_fit;
    ylo = min(flat_dB) - 0.5;
    yhi = max(flat_dB) + 0.5;

    if alphaInRange
        patch(ax, ...
            [alphaBand(1) alphaBand(2) alphaBand(2) alphaBand(1)], ...
            [ylo ylo yhi yhi], ...
            [0.88 0.93 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
            'HandleVisibility', 'off');
        hold(ax, 'on');
    end

    plot(ax, f_fit, flat_dB, 'b-', 'LineWidth', lw);
    yline(ax, 0, '--k', 'Aperiodic baseline', ...
        'LabelHorizontalAlignment', 'left', 'FontSize', fsTick - 1);

    if alphaInRange
        hold(ax, 'off');
        text(ax, mean(alphaBand), yhi * 0.85, '\alpha', ...
            'HorizontalAlignment', 'center', 'FontSize', fsTick + 1, ...
            'Color', [0.2 0.4 0.8], 'Interpreter', 'tex');
    end

    xlabel(ax, 'Frequency (Hz)',     'FontSize', fsAxis);
    ylabel(ax, 'Residual (dB)',      'FontSize', fsAxis);
    title(ax, 'Flattened Spectrum (Periodic Component)', ...
        'Interpreter', 'none', 'FontSize', fsAxis);
    xlim(ax, [fmin_fit fmax_fit]);
    ylim(ax, [ylo yhi]);
    set(ax, 'FontSize', fsTick); grid(ax, 'on');

    % [2,3] Per-trial fit RMSE time series
    ax = nexttile(tl, 6);
    plot(ax, tt, trialRMSE, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms, 'Color', [0.84 0.15 0.16]);
    if nT >= 5
        hold(ax, 'on');
        sm = smoothdata(trialRMSE, 'movmean', max(3, round(nT / 10)));
        plot(ax, tt, sm, '-', 'LineWidth', lw + 0.5, 'Color', [0.50 0.50 0.50]);
        legend(ax, {'Trial RMSE', 'Smoothed'}, 'Location', 'best', 'FontSize', fsTick - 1);
        hold(ax, 'off');
    end
    xlabel(ax, 'Trial',              'FontSize', fsAxis);
    ylabel(ax, 'Fit RMSE (dB)',      'FontSize', fsAxis);
    title(ax, 'Per-Trial Aperiodic Fit RMSE', 'Interpreter', 'none', 'FontSize', fsAxis);
    set(ax, 'FontSize', fsTick); grid(ax, 'on');

else
    % No PSD data supplied — informative placeholders for tiles 4-6
    placeholderMsg = 'PSD data not supplied (pass f and gaPxx)';
    titlesR2 = {'Mean PSD + Aperiodic Fit', 'Flattened Spectrum', 'Per-Trial Fit RMSE'};
    for k = 1:3
        ax = nexttile(tl, 3 + k);
        text(ax, 0.5, 0.5, placeholderMsg, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5], 'FontSize', fsTick);
        axis(ax, 'off');
        title(ax, titlesR2{k}, 'Interpreter', 'none', 'FontSize', fsAxis);
    end
end

exportgraphics(h, outPath, 'Resolution', dpi);
close(h);
end