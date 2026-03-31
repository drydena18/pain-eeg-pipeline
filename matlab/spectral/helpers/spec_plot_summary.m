function spec_plot_summary(outPath, f, gaPxx, featGA, fooofOut, alpha, cfg, subjid)
% SPEC_PLOT_SUMMARY  One-page per-subject summary figure
% V 1.3.0
%
% V1.3.0: calls both spec_plot_fooof_summary and spec_plot_fooof_aperiodic,
%         passing f and gaPxx to the latter so spectral panels are populated.
%
% V1.2.0: passed f and gaPxx to spec_plot_fooof_summary (now reverted;
%         spectral panels live in spec_plot_fooof_aperiodic instead).
%
% V1.1.0: panel 3 uses pre-stim interaction metrics with fallback to
%         full-epoch sf_balance / sf_logratio.

nTr = size(gaPxx, 2);

meanPSD = mean(gaPxx, 2, 'omitnan');

h = figure('Visible', 'off');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% ---- 1) Trial-averaged GA PSD ----
nexttile;
plot(f, 10*log10(meanPSD));
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title(sprintf('sub-%03d: GA PSD (trial-avg)', subjid), 'Interpreter', 'none');
xline(alpha.alpha_hz(1)); xline(alpha.alpha_hz(2));
xline(alpha.slow_hz(2), '--');
grid on;

% ---- 2) PAF CoG over trials ----
nexttile;
plot(1:nTr, featGA.paf_cog_hz, '-o', 'MarkerSize', 3);
yline(alpha.slow_hz(2), '--', '10 Hz', 'LabelHorizontalAlignment', 'left');
xlabel('Trial'); ylabel('Frequency (Hz)');
title('PAF (CoG) over trials', 'Interpreter', 'none');
grid on;

% ---- 3) Pre-stim interaction metrics (or fallback to full-epoch) ----
nexttile;
hasPreStim = isfield(featGA, 'bi_pre') && ~all(isnan(featGA.bi_pre));

if hasPreStim
    yyaxis left;
    plot(1:nTr, featGA.bi_pre, '-o', 'MarkerSize', 3); hold on;
    yline(0, ':', 'Color', [0.5 0.5 0.5]);
    ylabel('BI_{pre}');

    yyaxis right;
    plot(1:nTr, featGA.delta_erd, '-s', 'MarkerSize', 3);
    yline(0, ':', 'Color', [0.5 0.5 0.5]);
    ylabel('\DeltaERD');

    legend({'BI_{pre}', '\DeltaERD'}, 'Interpreter', 'tex', 'Location', 'best');
    title('Pre-stim BI & \DeltaERD over trials', 'Interpreter', 'tex');
else
    plot(1:nTr, featGA.sf_balance,  '-o', 'MarkerSize', 3); hold on;
    plot(1:nTr, featGA.sf_logratio, '-s', 'MarkerSize', 3);
    legend({'sf\_balance', 'sf\_logratio'}, 'Interpreter', 'none', 'Location', 'best');
    title('Alpha interaction over trials', 'Interpreter', 'none');
end
xlabel('Trial');
grid on;

% ---- 4) Placeholder tile (FOOOF saved to separate files) ----
nexttile;
axis off;
text(0.5, 0.5, {'FOOOF alpha peaks -> _fooof.png', 'FOOOF aperiodic  -> _fooof_aperiodic.png'}, ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'Color', [0.5 0.5 0.5]);

exportgraphics(h, outPath, 'Resolution', 200);
close(h);

% ---- FOOOF: alpha-peak / fit-quality summary (3x2, no PSD panels) ----
outFooofFig = strrep(outPath, '_summary.png', '_fooof.png');
spec_plot_fooof_summary(outFooofFig, fooofOut, subjid, cfg);

% ---- FOOOF: dedicated aperiodic figure (2x3, with PSD overlay) ----
outApFig = strrep(outPath, '_summary.png', '_fooof_aperiodic.png');
spec_plot_fooof_aperiodic(outApFig, fooofOut, subjid, cfg, f, gaPxx);
end