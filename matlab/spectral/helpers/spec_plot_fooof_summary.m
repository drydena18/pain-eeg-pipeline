function spec_plot_fooof_summary(outPath, fooofOut, subjid, cfg)

if nargin < 4
    cfg = struct();
end

% Defaults
dpi = 250;
fsTitle = 12;
fsAxis = 10;
fsTick = 9;
lw = 1.25;
mk = 'o';
ms = 3.5;

try
    if isfield(cfg, 'spectral') && isfield(cfg.spectral, 'qc') && isfield(cfg.spectral.qc, 'fooof_plot')
        fp = cfg.spectral.qc.fooof_plot;
        if isfield(fp, 'dpi'), dpi = fp.dpi; end
        if isfield(fp, 'fsTitle'), fsTitle = fp.fsTitle; end
        if isfield(fp, 'fsAxis'), fsAxis = fp.fsAxis; end
        if isfield(fp, 'fsTick'), fsTick = fp.fsTick; end
        if isfield(fp, 'lineWidth'), lw = fp.lineWidth; end
        if isfield(fp, 'markerSize'), ms = fp.markerSize; end
    end
catch
end

% Validate fooofOut
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

% Pull arrays
expnt = nan(nT, 1); offs = nan(nT, 1); r2 =  nan(nT, 1); err = nan(nT, 1);
acf = nan(nT, 1); apw = nan(nT, 1); abw = nan(nT, 1);

for t = 1:nT
    if isfield(tr(t), 'aperiodic_exponent'), expnt(t) = tr(t).aperiodic_exponent; end
    if isfield(tr(t), 'aperiodic_offset'),   offs(t)  = tr(t).aperiodic_offset; end
    if isfield(tr(t), 'r2'),                 r2(t)    = tr(t).r2; end
    if isfield(tr(t), 'error'),              err(t)   = tr(t).error; end
    if isfield(tr(t), 'alpha_cf'),           acf(t)   = tr(t).alpha_cf; end
    if isfield(tr(t), 'alpha_pw'),           apw(t)   = tr(t).alpha_pw; end
    if isfield(tr(t), 'alpha_bw'),           abw(t)   = tr(t).alpha_bw; end
end

tt = 1:nT;

% Layout
tl = tiledlayout(h, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Title
mainTitle = sprintf('sub-%03d FOOOF summary (trial-wise)', subjid);
if isfield(fooofOut, 'summary')
    s = fooofOut.summary;
    if isfield(s, 'fmin_hz') && isfield(s, 'fmax_hz')
        mainTitle = sprintf('%s | fit: %.1f - %.1f Hz', mainTitle, s.fmin_hz, s.fmax_hz);
    end
    if isfield(s, 'aperiodic_mode')
        mainTitle = sprintf('%s | %s', mainTitle, string(s.aperiodic_mode));
    end
end
title(tl, mainTitle, 'Interpreter', 'none', 'FontSize', fsTitle);

% 1) Aperiodic component
nexttile;
plot(tt, expnt, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Exponent', 'FontSize', fsAxis);
title('Aperiodic Exponent', 'Interpreter', 'none', 'FontSize', fsAxis);
set(gca, 'FontSize', fsTick); grid on;

% 2) Aperiodic offset
nexttile;
plot(tt, offs, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Offset', 'FontSize', fsAxis);
title('Aperiodic Offset', 'Interpreter', 'none', 'FontSize', fsTitle);
set(gca, 'FontSize', fsTick); grid on;

% 3) Alpha peak CF
nexttile;
plot(tt, acf, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Alpha CF (Hz)', 'FontSize', fsAxis);
title('Alpha Peak Center Frequency', 'Interpreter', 'none', 'FontSize', fsTitle);
set(gca, 'FontSize', fsTick); grid on;

% 4) Alpha peak power
nexttile;
plot(tt, apw, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Alpha peak power', 'FontSize', fsAxis);
title('Alpha Peak Power', 'Interpreter', 'none', 'FontSize', fsTitle);
set(gca, 'FontSize', fsTick); grid on;

% 5) Alpha peak bandwidth
nexttile;
plot(tt, abw, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Alpha BW (Hz)', 'FontSize', fsAxis);
title('Alpha peak bandwidth', 'Interpreter', 'none', 'FontSize', fsTitle);
set(gca, 'FontSize', fsTick); grid on;

% 6) Fit Quality
nexttile;
plot(tt, r2, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms); hold on;
plot(tt, err, ['-' mk], 'LineWidth', lw, 'MarkerSize', ms);
xlabel('Trial', 'FontSize', fsAxis); ylabel('Value', 'FontSize', fsAxis);
title('Fit quality', 'Interpreter', 'none', 'FontSize', fsTitle);
legend({'R^2', 'Error'}, 'Interpreter', 'none', 'Location', 'best');
set(gca, 'FontSize', fsTick); grid on;

exportgraphics(h, outPath, 'Resolution', dpi);
close(h);

end