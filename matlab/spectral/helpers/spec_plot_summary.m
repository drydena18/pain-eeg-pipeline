function spec_plot_summary(outPath, f, gaPxx, featGA, fooofOut, alpha, cfg, subjid)
% SPEC_PLOT_SUMMARY  One-page per-subject summary figure
% V 1.2.0
%
% V1.2.0: passes f and gaPxx to spec_plot_fooof_summary so that the
%         aperiodic fit overlay panels (7 and 8) are populated.
%
% V1.1.0: panel 3 updated to pre-stim interaction metrics (bi_pre, delta_erd)
%         with fallback to full-epoch sf_balance / sf_logratio.

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
    % Fallback: full-epoch metrics from V1
    plot(1:nTr, featGA.sf_balance,  '-o', 'MarkerSize', 3); hold on;
    plot(1:nTr, featGA.sf_logratio, '-s', 'MarkerSize', 3);
    legend({'sf\_balance', 'sf\_logratio'}, 'Interpreter', 'none', 'Location', 'best');
    title('Alpha interaction over trials', 'Interpreter', 'none');
end

xlabel('Trial');
grid on;

% ---- 4) placeholder (FOOOF summary is saved to its own file) ----
nexttile;
axis off;
text(0.5, 0.5, 'FOOOF → see _fooof.png', ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'Color', [0.5 0.5 0.5]);

exportgraphics(h, outPath, 'Resolution', 200);
close(h);

% ---- FOOOF detail figure (separate file, now with PSD overlay) ----
outFooofFig = strrep(outPath, '_summary.png', '_fooof.png');
spec_plot_fooof_summary(outFooofFig, fooofOut, subjid, cfg, f, gaPxx);
end