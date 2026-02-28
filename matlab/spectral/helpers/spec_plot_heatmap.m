function spec_plot_heatmap(outPath, X, chanLabels, ttl)
% X [chan x trial]
h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 2200 1400]);

imagesc(X);
axis tight;
set(gca, 'YDir', 'normal');

xlabel('Trial'); ylabel('Channel');
title(ttl, 'Interpreter', 'none');

nChan = numel(chanLabels);

step = 1;
if nChan > 40, step = 2; end
if nChan > 80, step = 4; end

yt = 1:step:nChan;
set(gca, 'YTick', yt, 'YTickLabel', chanLabels(yt));
set(gca, 'FontSize', 18);

colorbar;

exportgraphics(h, outPath, 'Resolution', 300);
close(h);
end