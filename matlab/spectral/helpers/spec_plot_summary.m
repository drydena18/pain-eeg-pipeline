function spec_plot_summary(outPath, f, gaPxx, featGA, fooofOut, alpha, cfg, subjid)
% SPEC_PLOT_SUMMARY  One-page per-subject summary figure
% V 1.1.0 — panel 3 updated to pre-stim interaction metrics (bi_pre, lr_pre,
%            cog_pre, delta_erd) with fallback to full-epoch sf_balance /
%            sf_logratio when pre-stim metrics are absent.

nTr = size(gaPxx, 2);

meanPSD = mean(gaPxx, 2, 'omitnan');

h = figure('Visible', 'off');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% ---- 1) Trial-avg GA PSD ----
nexttile;
plot(f, 10*log10(meanPSD));
xlabel('Frequency (Hz)'); ylabel('Power (dB)');
title(sprintf('sub-%03d: GA PSD (trial-avg)', subjid), 'Interpreter', 'none');
xline(alpha.alpha_hz(1)); xline(alpha.alpha_hz(2));
xline(alpha.slow_hz(2),  '--');
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
    plot(1:nTr, featGA.bi_pre,  '-o', 'MarkerSize', 3); hold on;
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

% ---- 4) FOOOF summary ----
outFooofFig = strrep(outPath, '_summary.png', '_fooof.png');
spec_plot_fooof_summary(outFooofFig, fooofOut, subjid, cfg);

exportgraphics(h, outPath, 'Resolution', 200);
close(h);
end