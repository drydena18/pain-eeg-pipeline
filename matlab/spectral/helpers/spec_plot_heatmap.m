function spec_plot_heatmap(outPath, X, chanLabels, ttl)
% X [chan x trial]
h = figure('Visible', 'off');
imagesc(X);
xlabel('Trial'); ylabel('Channel');
title(ttl, 'Interpreter', 'none');
set(gca, 'YTick', 1:numel(chanLabels), 'YTickLabel', chanLabels);
colorbar;
exportgraphics(h, outPath, 'Resolution', 200);
close(h);
end