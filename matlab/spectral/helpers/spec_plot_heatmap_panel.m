function spec_plot_heatmap_panel(outPath, featChan, chanLabels, subjid)

if isstring(chanLabels), chanLabels = cellstr(chanLabels); end

% Define the panel
items = {
    'paf_cog_hz', 'PAF CoG (Hz)';
    'pow_slow_alpha', 'Slow Alpha Power (8-10 Hz)';
    'pow_fast_alpha', 'Fast Alpha Power (10-12 Hz)';
    'pow_alpha_total', 'Total Alpha Power (8-12 Hz)';
    'rel_slow_alpha', 'Relative Slow Alpha Power';
    'rel_fast_alpha', 'Relative Fast Alpha Power';
    'sf_ratio', 'Slow/Fast Ratio';
    'sf_logratio', 'Slow/Fast Log Ratio';
    'sf_balance', 'Slow/Fast Balance';
    'slow_alpha_frac', 'Slow Alpha Fraction';
};

% Filter to existing fields
keep = false(size(items, 1), 1);
for k = 1:size(items, 1)
    fn = items{k, 1};
    if isfield(featChan, fn) && ~isempty(featChan.(fn))
        X = featChan.(fn);
        keep(k) = ismatrix(X) && ~isempty(X);
    end
end

items = items(keep, :);

if isempty(items)
    warning('spec_plot_heatmap_panel:NoFields', 'No valid featChan fields to plot.');
    return;
end

nPlots = size(items, 1);

if nPlots <= 4
    nR = 2; nC = 2;
elseif nPlots <= 6
    nR = 3; nC = 2;
elseif nPlots <= 9
    nR = 3; nC = 3;
else
    nR = 4; nC = 3;
end

nChan = numel(chanLabels);

% Figure sizing
figW = 2600;
figH = max(1400, 18*nChan);
h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 figW figH]);

tl = tiledlayout(h, nR, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('sub-%03d: Spectral Heatmaps (Channel x Trial)', subjid), 'Interpreter', 'none');

% Y-ticks
step = 1;
if nChan > 40, step = 2; end
if nChan > 80, step = 4; end
yt = 1:step:nChan;

for p = 1:nPlots
    fn = items{p, 1};
    ttl = items{p, 2};
    X = featChan.(fn);

    ax = nexttile(tl);
    imagesc(ax, X);
    axis(ax, 'tight');
    set(ax, 'YDir', 'normal');

    title(ax, ttl, 'Interpreter', 'none');

    xlabel(ax, 'Trial'); ylabel(ax, 'Channel');

    set(ax, 'YTick', yt, 'YTickLabel', chanLabels(yt));
    set(ax, 'FontSize', 11);

    cb = colorbar(ax);
    cb.Location = 'eastoutside';
end

exportgraphics(h, outPath, 'Resolution', 300);
close(h);
end