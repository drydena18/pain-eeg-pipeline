function spec_plot_heatmap_panel(outPath, featChan, chanLabels, subjid)
% SPEC_PLOT_HEATMAP_PANEL
% Channel x trial heatmaps for all spectral features
% V 1.1.0 - added pre-stim interaction metric fields


if isstring(chanLabels), chanLabels = cellstr(chanLabels); end

% Each row: {field_name, display_title}
% Fields absent from featChan are silently skipped
items = {
    % Full-epoch features (V1)
    'paf_cog_hz', 'PAF CoG (Hz)';
    'pow_slow_alpha', 'Slow Alpha Power (8–10 Hz)';
    'pow_fast_alpha', 'Fast Alpha Power (10–12 Hz)';
    'pow_alpha_total', 'Total Alpha Power (8–12 Hz)';
    'rel_slow_alpha', 'Relative Slow Alpha Power';
    'rel_fast_alpha', 'Relative Fast Alpha Power';
    'sf_ratio', 'Slow/Fast Ratio';
    'sf_logratio', 'Slow/Fast Log Ratio';
    'sf_balance', 'Slow-Fast Balance';
    'slow_alpha_frac', 'Slow Alpha Fraction';
    % Pre-stim interaction metrics (V2)
    'bi_pre', 'BI_{pre}: Pre-stim Balance Index';
    'lr_pre', 'LR_{pre}: Pre-stim Log Ratio';
    'cog_pre', 'CoG_{pre}: Pre-stim Center of Gravity (Hz)';
    'psi_cog', '\Psi_{CoG}: BI_{pre} x (CoG_{pre} - 10)';
    'erd_slow', 'ERD_{slow}: Slow-alpha ERD',
    'erd_fast', 'ERD_{fast}: Fast-alpha ERD';
    'delta_erd', '\DeltaERD: ERD_{slow} - ERD_{fast}';
    'p5_flag', 'p5 Flag (unstable pre-stim power)';
    % Phase (V2)
    'phase_slow_rad', 'Slow-alpha Phase at t = 0 (rad)';
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

nC = 5;
nR = ceil(nPlots / nC);

% Figure sizing (narrower + taller) 
figW = 900;                            % reduced width
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

    % Use a diverging colourmap for signed metrics centered on zero
    signed_fields = {'bi_pre', 'lr_pre', 'psi_cog', 'erd_slow', 'erd_fast', ...
        'delta_erd', 'sf_logratio', 'sf_balance', 'phase_slow_rad'};
    if ismember(fn, signed_fields)
        clim_abs = max(abs(X(:)), [], 'omitnan');
        if clim_abs > 0
            clim(ax, [-clim_abs, clim_abs]);
        end
        colormap(ax, spec_diverging_cmap());
    end

    colorbar(ax, 'eastoutside');
end

exportgraphics(h, outPath, 'Resolution',300);
close(h);

end

% ================================================================
% LOCAL: blue-white-red diverging colormap (64 levels)
% ================================================================
function cmap = spec_diverging_cmap()
n    = 32;
blue = [linspace(0.2, 1, n)', linspace(0.4, 1, n)', ones(n, 1)];
red  = [ones(n, 1), linspace(1, 0.2, n)', linspace(1, 0.2, n)'];
cmap = [blue; red];
end