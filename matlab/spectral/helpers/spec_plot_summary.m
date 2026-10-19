function spec_plot_summary(outPath, f, gaPxx, featGA, fooofOut, alpha, cfg, subjid)
% SPEC_PLOT_SUMMARY  One-page per-subject summary figure
% V 3.0.0
%
% V3.0.0 (feature set v3.1):
%   - delta_erd -> erd_asym (delta_ now always means post - pre). The old
%     name made the whole figure fail.
%   - Trial-wise FOOOF figures removed: FOOOF is now fitted once per
%     subject x channel on the trial-mean PSD (results in
%     sub-XXX_spectral_chan_summary.csv), so there is no per-trial series
%     to plot. fooofOut is accepted for call compatibility and ignored.
%   - Panel 4 now shows the per-trial pre-stimulus slow / fast alpha power
%     (GA) instead of a placeholder.
%
% Inputs:
%   gaPxx  : [nF x nTr] channel-average whole-window PSD per trial
%   featGA : struct of [1 x nTr] channel-average features

nTr = size(gaPxx, 2);
tr  = 1:nTr;

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
plot(tr, featGA.whole_paf_cog_hz, '-o', 'MarkerSize', 3);
yline(alpha.slow_hz(2), '--', '10 Hz', 'LabelHorizontalAlignment', 'left');
xlabel('Trial'); ylabel('Frequency (Hz)');
title('PAF (CoG) over trials', 'Interpreter', 'none');
grid on;

% ---- 3) Pre-stim balance + ERD asymmetry ----
nexttile;
hasPre = isfield(featGA, 'pre_sf_balance') && ~all(isnan(featGA.pre_sf_balance));
hasErd = isfield(featGA, 'erd_asym') && ~all(isnan(featGA.erd_asym));
if hasPre && hasErd
    yyaxis left;
    plot(tr, featGA.pre_sf_balance, '-o', 'MarkerSize', 3); hold on;
    yline(0, ':', 'Color', [0.5 0.5 0.5]);
    ylabel('BI_{pre}');

    yyaxis right;
    plot(tr, featGA.erd_asym, '-s', 'MarkerSize', 3);
    yline(0, ':', 'Color', [0.5 0.5 0.5]);
    ylabel('ERD_{asym} (slow - fast)');

    legend({'BI_{pre}', 'ERD_{asym}'}, 'Interpreter', 'tex', 'Location', 'best');
    title('Pre-stim BI & ERD asymmetry over trials', 'Interpreter', 'tex');
else
    plot(tr, featGA.whole_sf_balance,  '-o', 'MarkerSize', 3); hold on;
    plot(tr, featGA.whole_sf_logratio, '-s', 'MarkerSize', 3);
    legend({'sf\_balance', 'sf\_logratio'}, 'Interpreter', 'none', 'Location', 'best');
    title('Alpha interaction over trials', 'Interpreter', 'none');
end
xlabel('Trial');
grid on;

% ---- 4) Pre-stim slow / fast power over trials ----
nexttile;
if isfield(featGA, 'pre_pow_slow_alpha') && isfield(featGA, 'pre_pow_fast_alpha')
    plot(tr, featGA.pre_pow_slow_alpha, '-o', 'MarkerSize', 3); hold on;
    plot(tr, featGA.pre_pow_fast_alpha, '-s', 'MarkerSize', 3);
    legend({'slow (8-10 Hz)', 'fast (10-12 Hz)'}, 'Interpreter', 'none', 'Location', 'best');
    xlabel('Trial'); ylabel('Power');
    title('Pre-stim slow / fast alpha power', 'Interpreter', 'none');
    grid on;
else
    axis off;
end

exportgraphics(h, outPath, 'Resolution', 200);
close(h);
end