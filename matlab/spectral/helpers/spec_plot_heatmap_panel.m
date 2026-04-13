function spec_plot_heatmap_panel(outDir, featChan, chanLabels, subjid)
% SPEC_PLOT_HEATMAP_PANEL Channel x trial heatmaps split into three figures
% V 2.0.0
%
% V2.0.0 redesign:
%   - Three separare output tiles instead of one giant figure:
%       sub-XXX_heatmap_prestim.png     - full-epoch + pre-stim metrics
%       sub-XXX_heatmap_poststim.png    - ERD (post/pre metrics)
%       sub-XXX_heatmap_phase.png       - slow-alpha phase heatmap + rose plot
%   - Horizontal layout: panels tile left-to-right across 2 rows rather
%     than stacking vertically (was 2 wide x 7+ tall; now wide + 2 tall)
%   - Phase figure fixed: plarhistogram requires a PolarAxes parent.
%     The phase heatmap is shown left, the GA rose plot right.
%   - Fitst argument is now outDir (the figures/ folder) rather than a
%     single outPath. Filenames are derived automatically.
%
% Inputs:
%   outDir      : path to the subject figures directory
%   featChan    : struct of [nChan x nTr] spectral feature arrays
%   chanLabels  : cell or string array of channel labels
%   subjid      : integer subject ID

if isstring(chanLabels), chanLabels = cellstr(chanLabels); end

nChan = numel(chanLabels);
step = 1;
if nChan > 40, step = 2; end
if nChan > 80, step = 4; end
yt = 1:step:nChan;

% ---------------------------------------------------------------
% Field group definitions (field_name, display_title, signed?)
% Fields absent from featChan are silently skipped in each group
% ---------------------------------------------------------------

% --- Group 1: Pre-stim (full-epoch + pre-stimulus metrics) ---
preStimItems = {
    'paf_cog_hz', 'PAF CoG (Hz)', false;
    'pow_slow_alpha', 'Slow \alpha Power (8–10 Hz)', false;
    'pow_fast_alpha', 'Fast \alpha Power (10–12 Hz)', false;
    'pow_alpha_total', 'Total \alpha Power (8–12 Hz)', false;
    'rel_slow_alpha', 'Relative Slow \alpha', false;
    'rel_fast_alpha', 'Relative Fast \alpha', false;
    'sf_ratio', 'Slow/Fast Ratio', false;
    'sf_logratio', 'Slow/Fast Log Ratio', true;
    'sf_balance', 'Slow-Fast Balance', true;
    'slow_alpha_frac', 'Slow \alpha Fraction', false;
    'bi_pre', 'BI_{pre} | Pre-stim Balance Index', true;
    'lr_pre', 'LR_{pre} | Pre-stim Log Ratio', true;
    'cog_pre', 'CoG_{pre} | Pre-stim CoG (Hz)', false;
    'psi_cog', '\Psi_{cog} | BI x (CoG - 10)', true;
    'p5_flag', 'p5 Flag (unstable pre-stim power)', false;
};

% --- Group 2: Post-stim / ERD metrics ---
postStimItems = {
    'erd_slow', 'ERD_{slow} | Slow-\alpha ERD', true;
    'erd_fast', 'ERD_{fast} | Fast-\alpha ERD', true;
    'delta_erd', '\DeltaERD | ERD_{slow} - ERD_{fast}', true;
};

% ---------------------------------------------------------------
% Helper; filter item list to fields actually present in featChan
% ---------------------------------------------------------------
preStimItems = filter_items(preStimItems, featChan);
postStimItems = filter_items(postStimItems, featChan);
hasPhase = isfield(featChan, 'phase_slow_rad') && ...
    ~isempty(featChan.phase_slow_rad);

% ================================================================
% FIGURE 1: Pre-stim heatmpas (horizontal layout)
% ================================================================
if ~isempty(preStimItems)
    nP = size(preStimItems, 1);
    nC = max(1, ceil(nP / 2));
    nR = ceil(nP / nC);

    % Wide figure: each panel ~320px wide; height from channel count
    panelH = max(250, 14*nChan + 80);
    figW = nC * 750;
    figH = nR * panelH;

    h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [50 50 figW figH]);
    tl = tiledlayout(h, nR, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('sub-%03d Pre-stim Spectral Heatmaps [chan x trial]', subjid), ...
        'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');

    for p = 1:nP
        fn = preStimItems{p, 1};
        ttl = preStimItems{p, 2};
        signed = preStimItems{p, 3};
        X = featChan.(fn);

        ax = nexttile(tl);
        imagesc(ax, X);
        axis(ax, 'tight');
        set(ax, 'YDir', 'normal');
        set(ax, 'YTick', yt, 'YTickLabel', chanLabels(yt), 'FontSize', 8);
        title(ax, ttl, 'Interpreter', 'tex', 'FontSize', 9);
        xlabel(ax, 'Trial', 'FontSize', 8);
        ylabel(ax, 'Channel', 'FOntSize', 8);

        if signed
            clim_abs = max(abs(X(:)), [], 'omitnan');
            if clim_abs > 0, clim(ax, [-clim_abs clim_abs]); end
            colormap(ax, spec_diverging_cmap());
        else
            colormap(ax, parula);
        end
        colorbar(ax, 'Location', 'eastoutside');
    end

    ourPre = fullfile(outDir, sprintf('sub-%03d_heatmap_prestim.png', subjid));
    exportgraphics(h, ourPre, 'Resolution', 200);
    close(h);
end

% ================================================================
% FIGURE 2: Post-stim / ERD Heatmaps (single row)
% ================================================================
if ~isempty(postStimItems)
    nP = size(postStimItems, 1);
    nC = nP;
    nR = 1;

    panelH = max(350, 14*nChan + 80);
    figW = nC * 750;
    figH = panelH;

    h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [50 50 figW figH]);
    tl = tiledlayout(h, nR, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('sub-%03d Post-stim ERD Heatmaps [chan x trial]', subjid), ...
        'Interpreter', 'tex', 'FontSize', 11, 'FontWeight', 'bold');
    
    for p = 1:nP
        fn = postStimItems{p, 1};
        ttl = postStimItems{p, 2};
        signed = postStimItems{p, 3};
        X = featChan.(fn);

        ax = nexttile(tl);
        imagesc(ax, X);
        axis(ax, 'tight');
        set(ax, 'YDir', 'normal');
        set(ax, 'YTick', yt, 'YTickLabel', chanLabels(yt), 'FontSize', 8);
        title(ax, ttl, 'Interpreter', 'tex', 'FontSize', 9);
        xlabel(ax, 'Trial', 'FontSize', 8);
        ylabel(ax, 'Channel', 'FontSize', 8);

        if signed
            clim_abs = max(abs(X(:)), [], 'omitnan');
            if clim_abs > 0, clim(ax, [-clim_abs clim_abs]); end
            colormap(ax, spec_diverging_cmap());
        else
            colormap(ax, parula);
        end
        colormap(ax, 'Location','eastoutside');
    end

    outPost = fullfile(outDir, sprintf('sub-%03d_heatmap_poststim.png', subjid));
    exportgraphics(h, outPost, 'Resolution', 200);
    close(h);
end

% ================================================================
% FIGURE 3: Phase heatmap + GA rose plot
%
% FIX: plorhistogram requires a PolarAxes parent. nexttile returns
% a regular Axes, so we capture its Position, delete it, and create a
% polaraxes in the same screen rectangle.
% ================================================================
if hasPhase
    phi_mat = featChan.phase_slow_rad;

    % GA phase = circular mean across channels per trial
    phi_ga = angle(mean(exp(1i * phi_mat), 1, 'omitnan'));
    phi_ga_valid = phi_ga(~isnan(phi_ga));

    panelH = max(400, 14*nChan + 80);
    figW = 900;
    figH = panelH;

    h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [50 50 figW figH]);
    tl = tiledlayout(h, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('sub-%03d Slow-alpha Phase at t = 0', subjid), ...
        'Interpreter', 'tex', 'FontSize', 11, 'FontWeight', 'bold');

    % --- Left: channel x trial heatmap ---
    ax1 = nexttile(tl, 1);
    imagesc(ax1, phi_mat);
    axis(ax1, phi_mat);
    set(ax1, 'YDir', 'normal');
    set(ax1, 'YTick', yt, 'YTickLabel', chanLabels(yt), 'FontSize', 9);
    xlabel(ax1, 'Trial', 'FontSize', 8);
    ylabel(ax1, 'Channel', 'FontSize', 8);
    clim_abs = pi;
    clim(ax1, [-clim_abs clim_abs]);
    colormap(ax1, hsv(256));
    cb = colorbar(ax1, 'Location','eastoutside');
    cb.Ticks = [-pi -pi/2 0 pi/2 pi];
    cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};

    % --- Right: GA rose plot (polar histogram) ---
    % nexttile gives a regular Axes; replace with polaraxes in same position
    ax2_dummy = nexttile(tl, 2);
    pos2 = ax2_dummy.Position;
    delete(ax2_dummy);
    ax2 = polaraxes(h, 'Position', pos2);

    if ~isempty(phi_ga_valid)
        polarhistogram(ax2, phi_ga_valid, 24, ...
            'Normalization', 'probability', ...
            'FaceColor', [0.13 0.47 0.71], 'FaceAlpha', 0.65);
        ax2.ThetaZeroLocation = 'top';
        ax2.ThetaDir = 'clockwise';
        title(ax2, 'GA Phase Distribution (rose)', 'Interpreter', 'none', 'FontSize', 9);
    else
        title(ax2, 'GA Phase: no valid data', 'Interpreter', 'none', 'FontSize', 9);
    end

    outPh = fullfile(outDir, sprintf('sub-%03d_heatmap_phase.png', subjid));
    exportgraphics(h, outPh, 'Resolution', 200);
    close(h);
end

end

% ================================================================
% Local: keep only item rows where field existsin featChan
% ================================================================
function items = filter_items(items, featChan)
keep = false(size(items, 1), 1);
for k = 1:size(items, 1)
    fn = items{k, 1};
    if isfield(featChan, fn) && ~isempty(featChan.(fn)) && ismatrix(featChan.(fn))
        keep(k) = true;
    end
end
items = items(keep, :);
end

% ================================================================
% Local: blue-white-red diverging colormap (64 levels)
% ================================================================
function cmap = spec_diverging_cmap()
n = 32;
blue = [linspace(0.2, 1, n)', linspace(0.4, 1, n)', ones(n, 1)];
red = [ones(n, 1), linspace(1, 0.2, n)', linspace(1, 0.2, n)'];
cmap = [blue; red];
end