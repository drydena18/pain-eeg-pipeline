function spec_plot_heatmap_panel(outPath, featChan, chanLabels, subjid)

if isstring(chanLabels), chanLabels = cellstr(chanLabels); end

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

% --- Keep only valid fields ---
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
nChan  = numel(chanLabels);

% Layout strategy: max 2 columns
if nPlots <= 2
    nC = 1;
else
    nC = 2;
end
nR = ceil(nPlots / nC);

% Figure sizing (narrower + taller)
figW = 1100;                           % reduced width
figH = max(1200, 30*nChan + 250*nR);   % taller scaling

h = figure('Visible','off','Units','pixels','Position',[100 100 figW figH]);

tl = tiledlayout(h, nR, nC, ...
    'TileSpacing','compact', ...
    'Padding','compact');

title(tl, sprintf('sub-%03d: Spectral Heatmaps (Channel x Trial)', subjid), ...
    'Interpreter','none');

% Y tick scaling
step = 1;
if nChan > 40, step = 2; end
if nChan > 80, step = 4; end
yt = 1:step:nChan;

for p = 1:nPlots
    fn  = items{p,1};
    ttl = items{p,2};
    X   = featChan.(fn);

    ax = nexttile(tl);
    imagesc(ax, X);
    axis(ax,'tight');
    set(ax,'YDir','normal');

    title(ax, ttl, 'Interpreter','none');
    xlabel(ax,'Trial');
    ylabel(ax,'Channel');

    set(ax,'YTick',yt,'YTickLabel',chanLabels(yt),'FontSize',11);

    colorbar(ax,'eastoutside');
end

exportgraphics(h, outPath, 'Resolution',300);
close(h);

end