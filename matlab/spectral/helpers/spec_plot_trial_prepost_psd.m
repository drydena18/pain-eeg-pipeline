function spec_plot_trial_prepost_psd(outPath, fPre, pPre, fPost, pPost, chanLabels, subjid, trialIdx, legendMax)
% Pre/post: pPre, pPost are [Chan x Freq]

if isstring(chanLabels), chanLabels = cellstr(chanLabels); end
chanLabels = chanLabels(:); % force column cell

h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 1800 800]);

tl = tiledlayout(h, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('sub-%03d trial-%03d: Pre vs Post PSD (channels overlaid)', subjid, trialIdx), 'Interpreter', 'none');

nChan = size(pPre, 1);
idxLeg = 1:min(nChan, legendMax);

% ---------- PRE ----------
ax1 = nexttile(tl);

% Plot only subset channels
plot(ax1, fPre, pPre(idxLeg, :)');
set(ax1, 'YScale', 'log');
grid(ax1, 'on');
xlabel(ax1, 'Frequency (Hz)');
ylabel(ax1, 'Power (linear, log y)');
title(ax1, 'Pre-stim PSD', 'Interpreter', 'none');
xlim(ax1, [min(fPre) max(fPre)]);

hold(ax1, 'on');
plot(ax1, fPre, median(pPre, 1, 'omitnan'), 'LineWidth', 2);
hold(ax1, 'off');

legNames = [chanLabels(idxLeg).', {'MEDIAN'}]; 
legend(ax1, legNames, 'Location', 'northeastoutside');

% ---------- POST ----------
ax2 = nexttile(tl);

plot(ax2, fPost, pPost(idxLeg, :)');
set(ax2, 'YScale', 'log');
grid(ax2, 'on');
xlabel(ax2, 'Frequency (Hz)');
ylabel(ax2, 'Power (linear, log y)');
title(ax2, 'Post-stim PSD', 'Interpreter', 'none');
xlim(ax2, [min(fPost) max(fPost)]);

hold(ax2, 'on');
plot(ax2, fPost, median(pPost, 1, 'omitnan'), 'LineWidth', 2);
hold(ax2, 'off');

exportgraphics(h, outPath, 'Resolution', 300);
close(h);

end