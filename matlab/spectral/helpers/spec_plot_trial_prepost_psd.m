function spec_plot_trial_prepost_psd(outPath, fPre, pPre, fPost, pPost, chanLabels, subjid, trialIdx, legendMax)
% Pre/post: [Chan x Freq]

if isstring(chanLabels), chanLabels = cellstr(chanLabels); end

h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 1600 800]);

tl = tiledlayout(h, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('sub-%03d trial-%03d: Pre vs Post PSD (channels overlaid)', subjid, trialIdx), 'Interpreter', 'none');

% Decide which channels to include
nChan = size(pPre, 1);
idxLeg = 1:min(nChan, legendMax);

% ---- PRE ----
ax1 = nexttile(tl);
plot(ax1, fPre, pPre');
set(ax1, 'YScale', 'log');
grid(ax1, 'on');
xlabel(ax1, 'Frequency (Hz)'); ylabel(ax1, 'Power log scale (dB)');
title(ax1, 'Pre-stim PSD', 'Interpreter', 'none');
xlim(ax1, [min(fPre) max(fPre)]);

% Median across all channels
hold(ax1, 'on');
plot(ax1, fPre, median(pPre, 1, 'omitnan'), 'LineWidth', 2);
hold(ax1, 'off');

% Legend only for subset + median
legNames = [chanLabels(idxLeg), {'MEDIAN'}];
legend(ax1, legNames{:}, 'Location', 'northeastoutside');

% ---- POST ----
ax2 = nexttile(tl);
plot(ax2, fPost, pPost');
set(ax2, 'YScale', 'log');
grid(ax2, 'on');
xlabel(ax2, 'Frequency (Hz)'); ylabel(ax2, 'Power log scale (dB)');
title(ax2, 'Post-stim PSD', 'Interpreter', 'none');
xlim(ax2, [min(fPost) max(fPost)]);

hold(ax2, 'on');
plot(ax2, fPost, median(pPost, 1, 'omitnan'), 'LineWidth', 2);
hold(ax2, 'off');

exportgraphics(h, outPath, 'Resolution', 300);
close(h);

end