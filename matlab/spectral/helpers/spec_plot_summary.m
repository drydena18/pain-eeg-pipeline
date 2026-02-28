function spec_plot_summary(outPath, f, gaPxx, featGA, fooofOut, alpha, cfg, subjid)
% One-page summary: GA PSD, PAF over trials, alpha interaction over trials, FOOOF example
nTr = size(gaPxx, 2);

% Trial averaged GA PSD
meanPSD = mean(gzPxx, 2, 'omitnan');

% Choose two alpha interaction traces to show
y1 = featGA.sf_balance;
y2 = featGA.sf_logratio;

h = figure('Visible', 'off');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% 1) Trial-avg GA PSD
nexttile;
plot(f, 10*log10(meanPSD));
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title(sprintf('sub-%03d: GA PSD (trial-avg)', subjid), 'Interpreter', 'none');
xline(alpha.alpha_hz(1)); xline(alpha.alpha_hz(2));
xline(alpha.slow_hz(1), '--'); xline(alpha.slow_hz(2), '--');
xline(alpha.fast_hz(1), '--'); xline(alpha.fast_hz(2, '--');

% 2) PAF over trials
nexttile;
plot(1:nTr, featGA.paf_cog_hz, '-o');
xlabel('Trial'); ylabel('Frequency (Hz)');
title('PAF (CoG) over trials', 'Interpreter', 'none');

% 3) Alpha interaction metrics
nexttile;
plot(1:nTr, y1, '-o'); hold on;
plot(1:nTr, y2, '-o');
xlabel('Trial'); ylabel('a.u.');
legend({'sf\balance', 'sf\logratio'}, 'Interpreter', 'none', 'Location', 'best');
title('Alpha interaction over trials', 'Interpreter', 'none');

% 4) FOOOF summary (if available)
nexttile;
axis off;
if ~isempty(fooofOut) && isfield(fooofOut, 'summary')
    s = fooofOut.summary;
    txt = sprintf(['FOOOF (GA per-trial)\n', ...
        'f range: %.1f-%.1f Hz\n', ...
        'aperiodic_mode: %s\n', ...
        'alpha band: %.1f-%.1f Hz\n', ...
        'trials fit: %d\n'], ...
        s.fmin_hz, s.fmax_hz, s.aperiodic_mode, s.alpha_band_hz(1), s.alpha_band_hz(2), s.n_trials);
    text(0, 1, txt, 'VerticalAlignment', 'top', 'Interpreter', 'none');
else
    text(0, 1, 'FOOOF: not available / failed', 'VerticalAlignment', 'top', 'Interpreter', 'none');
end

exportgraphics(h, outPath, 'Resolution', 200);
close(h);
end