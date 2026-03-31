function spec_plot_fooof_summary(outPath, fooofOut, subjid, cfg, f, gaPxx)
% SPEC_PLOT_FOOOF_SUMMARY  Trial-wise FOOOF parameter plots + aperiodic fit
% V 1.1.0
%
% V1.1.0 adds two panels (7 and 8) showing the aperiodic model fit
% overlaid on the raw power spectral density (PSD):
%
%   Panel 7 — Mean raw PSD (log10 scale, dB) across trials with the
%             trial-averaged aperiodic fit overlaid.  The fit covers only
%             the FOOOF frequency range [fmin, fmax].  This makes the 1/f
%             structure and its model accuracy immediately visible.
%
%   Panel 8 — Trial-averaged "flattened" PSD: raw PSD minus the aperiodic
%             fit in log10 space.  Residuals > 0 are peaks (periodic
%             components); values near 0 indicate the aperiodic-only
%             region.  The alpha-band window is shaded for reference.
%
% Aperiodic fit reconstruction (fixed mode):
%   log10(power_ap) = offset - exponent * log10(f)
%   power_ap        = 10 ^ (offset - exponent * log10(f))
%
% For knee mode (aperiodic_knee field present and non-NaN):
%   log10(power_ap) = offset - log10(knee + f ^ exponent)
%
% Layout: 4 rows x 2 columns (tiles 1-6 unchanged from V1.0; tiles 7-8 new)
%
% Inputs:
%   outPath  : full path for the saved PNG
%   fooofOut : struct from spec_run_fooof_python / spec_fill_fooof_alpha
%   subjid   : integer subject ID
%   cfg      : pipeline cfg struct (optional; for plot styling overrides)
%   f        : [1 x nFreq] frequency vector (Hz) — the full-range PSD freqs
%              Pass [] to skip panels 7 and 8.
%   gaPxx    : [nFreq x nTrials] GA PSD across channels (linear power, μV²/Hz)
%              Pass [] to skip panels 7 and 8.

if nargin < 4, cfg   = struct(); end
if nargin < 5, f     = [];       end
if nargin < 6, gaPxx = [];       end

% ---------------------------------------------------------------
% Styling defaults (overridable via cfg.spectral.qc.fooof_plot)
% ---------------------------------------------------------------
dpi     = 250;
fsTitle = 12;
fsAxis  = 10;
fsTick  = 9;
lw      = 1.25;
mk      = 'o';
ms      = 3.5;

try
    if isfield(cfg, 'spectral') && isfield(cfg.spectral, 'qc') && isfield(cfg.spectral.qc, 'fooof_plot')
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
% Validate FOOOF output
% ---------------------------------------------------------------
hasFooof = isstruct(fooofOut) && isfield(fooofOut, 'trials') && ~isempty(fooofOut.trials);

h = figure('Visible', 'off');
set(h, 'Color', 'w');

if ~hasFooof
    axis off;
    text(0, 0.8, sprintf('sub-%03d FOOOF: not available / failed', subjid), ...
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
acf   = nan(nT, 1);
apw   = nan(nT, 1);
abw   = nan(nT, 1);

for t = 1:nT
    if isfield(tr(t), 'aperiodic_exponent'), expnt(t) = tr(t).aperiodic_exponent; end
    if isfield(tr(t), 'aperiodic_offset'),   offs(t)  = tr(t).aperiodic_offset;   end
    if isfield(tr(t), 'aperiodic_knee'),     knee(t)  = tr(t).aperiodic_knee;     end
    if isfield(tr(t), 'r2'),                 r2(t)    = tr(t).r2;                 end
    if isfield(tr(t), 'error'),              err(t)   = tr(t).error;              end
    if isfield(tr(t), 'alpha_cf'),           acf(t)   = tr(t).alpha_cf;           end
    if isfield(tr(t), 'alpha_pw'),           apw(t)   = tr(t).alpha_pw;           end
    if isfield(tr(t), 'alpha_bw'),           abw(t)   = tr(t).alpha_bw;           end
end

tt = (1:nT)';

% ---------------------------------------------------------------
% Figure layout: 4 rows x 2 cols (panels 7-8 are new)
% ---------------------------------------------------------------
tl = tiledlayout(h, 4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

mainTitle = sprintf('sub-%03d FOOOF summary (trial-wise)', subjid);
if isfield(fooofOut, 'summary')
    s = fooofOut.summary;
    if isfield(s, 'fmin_hz') && isfield(s, 'fmax_hz')
        mainTitle = sprintf('%s | fit: %.1f–%.1f Hz', mainTitle, s.fmin_hz, s.fmax_hz);
    end
    if isfield(s, 'aperiodic_mode')
        mainTitle = sprintf('%s | %s', mainTitle, string(s.aperiodic_mode));
    end
end
title(tl, mainTitle, 'Interpreter', 'none', 'FontSize', fsTitle);

% ---------------------------------------------------------------
% Panels 1-6: scalar parameter time-series (unchanged from V1.0)
% ---------------------------------------------------------------
nexttile;
plot(tt, expnt, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Exponent', 'FontSize', fsAxis);
title('Aperiodic Exponent', 'Interpreter', 'none', 'FontSize', fsAxis);
set(gca, 'FontSize', fsTick); grid on;

nexttile;
plot(tt, offs, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Offset', 'FontSize', fsAxis);
title('Aperiodic Offset', 'Interpreter', 'none', 'FontSize', fsAxis);
set(gca, 'FontSize', fsTick); grid on;

nexttile;
plot(tt, acf, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Alpha CF (Hz)', 'FontSize', fsAxis);
title('Alpha Peak Center Frequency', 'Interpreter', 'none', 'FontSize', fsAxis);
set(gca, 'FontSize', fsTick); grid on;

nexttile;
plot(tt, apw, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Alpha peak power', 'FontSize', fsAxis);
title('Alpha Peak Power', 'Interpreter', 'none', 'FontSize', fsAxis);
set(gca, 'FontSize', fsTick); grid on;

nexttile;
plot(tt, abw, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Alpha BW (Hz)', 'FontSize', fsAxis);
title('Alpha Peak Bandwidth', 'Interpreter', 'none', 'FontSize', fsAxis);
set(gca, 'FontSize', fsTick); grid on;

nexttile;
plot(tt, r2,  ['-' mk], 'LineWidth', lw, 'MarkerSize', ms); hold on;
plot(tt, err, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Value', 'FontSize', fsAxis);
title('Fit Quality', 'Interpreter', 'none', 'FontSize', fsAxis);
legend({'R^2', 'Error'}, 'Interpreter', 'none', 'Location', 'best');
set(gca, 'FontSize', fsTick); grid on;

% ---------------------------------------------------------------
% Panels 7-8: aperiodic fit overlaid on mean PSD
% Requires f and gaPxx from the caller.  If absent or FOOOF range
% cannot be determined, show a placeholder.
% ---------------------------------------------------------------
hasPSD = ~isempty(f) && ~isempty(gaPxx) && size(gaPxx, 1) == numel(f);

fmin_fit = NaN; fmax_fit = NaN;
if isfield(fooofOut, 'summary')
    if isfield(fooofOut.summary, 'fmin_hz'), fmin_fit = fooofOut.summary.fmin_hz; end
    if isfield(fooofOut.summary, 'fmax_hz'), fmax_fit = fooofOut.summary.fmax_hz; end
end

if hasPSD && ~isnan(fmin_fit) && ~isnan(fmax_fit)
    % Mean raw PSD across trials (linear power)
    meanPsd = mean(gaPxx, 2, 'omitnan');    % [nFreq x 1]

    % Restrict to FOOOF fitting range
    idxFit   = (f >= fmin_fit) & (f <= fmax_fit);
    f_fit    = f(idxFit);
    psd_fit  = meanPsd(idxFit);

    % Compute mean aperiodic parameters (omit failed trials = NaN)
    mean_offs = mean(offs, 'omitnan');
    mean_exp  = mean(expnt, 'omitnan');
    mean_knee = mean(knee, 'omitnan');   % NaN for fixed mode

    % Reconstruct aperiodic fit
    %   fixed mode: log10(ap) = offset - exponent * log10(f)
    %   knee  mode: log10(ap) = offset - log10(knee + f^exponent)
    useKnee = ~isnan(mean_knee) && mean_knee > 0;

    log10f = log10(f_fit(:));
    if useKnee
        log10_ap = mean_offs - log10(mean_knee + f_fit(:) .^ mean_exp);
    else
        log10_ap = mean_offs - mean_exp .* log10f;
    end

    % Convert to dB (10*log10(power))
    psd_fit_dB = 10 * log10(max(psd_fit(:), 1e-30));
    ap_fit_dB  = 10 * log10(max(10 .^ log10_ap, 1e-30));

    % Flattened (periodic) component = raw - aperiodic in log10 space
    flat_dB = psd_fit_dB - ap_fit_dB;

    % ---- Panel 7: raw PSD + aperiodic overlay ----
    ax7 = nexttile;
    plot(ax7, f_fit, psd_fit_dB, 'k-',  'LineWidth', lw,      'DisplayName', 'Raw PSD');
    hold(ax7, 'on');
    plot(ax7, f_fit, ap_fit_dB,  'r--', 'LineWidth', lw + 0.5,'DisplayName', 'Aperiodic fit');
    hold(ax7, 'off');

    xlabel(ax7, 'Frequency (Hz)', 'FontSize', fsAxis);
    ylabel(ax7, 'Power (dB re μV²/Hz)', 'FontSize', fsAxis);
    modeStr = 'fixed';
    if useKnee, modeStr = 'knee'; end
    title(ax7, sprintf('Mean PSD + Aperiodic Fit (%s)', modeStr), ...
        'Interpreter', 'none', 'FontSize', fsAxis);
    legend(ax7, 'Location', 'northeast', 'FontSize', fsTick - 1);
    xlim(ax7, [fmin_fit fmax_fit]);
    set(ax7, 'FontSize', fsTick);
    grid(ax7, 'on');

    % ---- Panel 8: flattened PSD (periodic component) ----
    ax8 = nexttile;

    % Shade alpha band if it falls within the fit range
    alphaBand = [8 12];
    if isfield(fooofOut, 'summary') && isfield(fooofOut.summary, 'alpha_band_hz')
        ab = fooofOut.summary.alpha_band_hz;
        if numel(ab) == 2, alphaBand = ab; end
    end
    alphaInRange = alphaBand(1) >= fmin_fit && alphaBand(2) <= fmax_fit;

    if alphaInRange
        patch(ax8, ...
            [alphaBand(1) alphaBand(2) alphaBand(2) alphaBand(1)], ...
            [min(flat_dB)-1  min(flat_dB)-1  max(flat_dB)+1  max(flat_dB)+1], ...
            [0.90 0.95 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
            'HandleVisibility', 'off');
        hold(ax8, 'on');
    end

    plot(ax8, f_fit, flat_dB, 'b-', 'LineWidth', lw);
    yline(ax8, 0, '--k', 'Aperiodic', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
    if alphaInRange
        hold(ax8, 'off');
        text(ax8, mean(alphaBand), max(flat_dB)*0.85, 'alpha', ...
            'HorizontalAlignment', 'center', 'FontSize', fsTick - 1, 'Color', [0.2 0.4 0.8]);
    end

    xlabel(ax8, 'Frequency (Hz)', 'FontSize', fsAxis);
    ylabel(ax8, 'Residual power (dB)', 'FontSize', fsAxis);
    title(ax8, 'Flattened PSD (Periodic Component)', ...
        'Interpreter', 'none', 'FontSize', fsAxis);
    xlim(ax8, [fmin_fit fmax_fit]);
    set(ax8, 'FontSize', fsTick);
    grid(ax8, 'on');

else
    % No PSD data — show informative placeholders
    ax7 = nexttile;
    text(ax7, 0.5, 0.5, 'PSD data not supplied (pass f and gaPxx)', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Color', [0.5 0.5 0.5], 'FontSize', fsTick);
    axis(ax7, 'off');
    title(ax7, 'Mean PSD + Aperiodic Fit', 'Interpreter', 'none', 'FontSize', fsAxis);

    ax8 = nexttile;
    text(ax8, 0.5, 0.5, 'PSD data not supplied (pass f and gaPxx)', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Color', [0.5 0.5 0.5], 'FontSize', fsTick);
    axis(ax8, 'off');
    title(ax8, 'Flattened PSD (Periodic Component)', 'Interpreter', 'none', 'FontSize', fsAxis);
end

exportgraphics(h, outPath, 'Resolution', dpi);
close(h);
end