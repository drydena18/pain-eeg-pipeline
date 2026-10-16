function spec_plot_heatmap_panel(outDir, featChan, chanLabels, subjid)
% SPEC_PLOT_HEATMAP_PANEL Channel x trial heatmaps split into three figures
% V 2.2.0
%
% V 2.2.0 changes vs V 2.1.0 (feature set v3.1 names):
%   - delta_erd -> erd_asym (ERD_slow - ERD_fast).
%   - whole_slow_alpha_frac removed (dropped from the feature set: it equals
%     (1 + sf_balance) / 2 exactly).
%   - phase_slow_rad -> slow_phase. The phase figure was silently skipped
%     because the old field no longer existed.
%   Fields still absent from featChan are skipped, as before.
%
% V 2.1.0 changes vs V 2.0.0:
%   - Field renames to match the whole_/pre_/post_/delta_ metric scheme.
%   - Added erd_pow_alpha_total to Group 2.
%
% Outputs (in outDir):
%   sub-XXX_heatmap_prestim.png     - whole-epoch + pre-stim metrics
%   sub-XXX_heatmap_poststim.png    - ERD metrics
%   sub-XXX_heatmap_phase.png       - slow-alpha phase heatmap + rose plot
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
% ---------------------------------------------------------------

% --- Group 1: whole-epoch + pre-stimulus metrics ---
preStimItems = {
    'whole_paf_cog_hz', 'PAF CoG (Hz)', false;
    'whole_pow_slow_alpha', 'Slow \alpha Power (8-10 Hz)', false;
    'whole_pow_fast_alpha', 'Fast \alpha Power (10-12 Hz)', false;
    'whole_pow_alpha_total', 'Total \alpha Power (8-12 Hz)', false;
    'whole_rel_slow_alpha', 'Relative Slow \alpha', false;
    'whole_rel_fast_alpha', 'Relative Fast \alpha', false;
    'whole_sf_ratio', 'Slow/Fast Ratio', false;
    'whole_sf_logratio', 'Slow/Fast Log Ratio', true;
    'whole_sf_balance', 'Slow-Fast Balance', true;
    'pre_sf_balance', 'BI_{pre} | Pre-stim Balance Index', true;
    'pre_sf_logratio', 'LR_{pre} | Pre-stim Log Ratio', true;
    'pre_paf_cog_hz', 'CoG_{pre} | Pre-stim CoG (Hz)', false;
    'pre_psi_cog', '\Psi_{cog} | BI x (CoG - 10)', true;
    'p5_flag', 'p5 Flag (unstable pre-stim power)', false;
};

% --- Group 2: Post-stim / ERD metrics ---
postStimItems = {
    'erd_slow', 'ERD_{slow} | Slow-\alpha ERD', true;
    'erd_fast', 'ERD_{fast} | Fast-\alpha ERD', true;
    'erd_pow_alpha_total', 'ERD_{total} | Total-\alpha ERD', true;
    'erd_asym', 'ERD_{asym} | ERD_{slow} - ERD_{fast}', true;
};

preStimItems = filter_items(preStimItems, featChan);
postStimItems = filter_items(postStimItems, featChan);
hasPhase = isfield(featChan, 'slow_phase') && ~isempty(featChan.slow_phase);

% ================================================================
% FIGURE 1: Pre-stim heatmaps (horizontal layout)
% ================================================================
if ~isempty(preStimItems)
    nP = size(preStimItems, 1);
    nC = max(1, ceil(nP / 2));
    nR = ceil(nP / nC);

    panelH = max(250, 14*nChan + 80);
    figW = nC * 750;
    figH = nR * panelH;

    h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [50 50 figW figH]);
    tl = tiledlayout(h, nR, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('sub-%03d Pre-stim Spectral Heatmaps [chan x trial]', subjid), ...
        'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');

    for p = 1:nP
        draw_heatmap(nexttile(tl), featChan.(preStimItems{p, 1}), preStimItems{p, 2}, ...
            preStimItems{p, 3}, chanLabels, yt);
    end

    exportgraphics(h, fullfile(outDir, sprintf('sub-%03d_heatmap_prestim.png', subjid)), 'Resolution', 200);
    close(h);
end

% ================================================================
% FIGURE 2: Post-stim / ERD heatmaps (single row)
% ================================================================
if ~isempty(postStimItems)
    nP = size(postStimItems, 1);

    panelH = max(350, 14*nChan + 80);
    h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [50 50 nP*750 panelH]);
    tl = tiledlayout(h, 1, nP, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('sub-%03d Post-stim ERD Heatmaps [chan x trial]', subjid), ...
        'Interpreter', 'tex', 'FontSize', 11, 'FontWeight', 'bold');

    for p = 1:nP
        draw_heatmap(nexttile(tl), featChan.(postStimItems{p, 1}), postStimItems{p, 2}, ...
            postStimItems{p, 3}, chanLabels, yt);
    end

    exportgraphics(h, fullfile(outDir, sprintf('sub-%03d_heatmap_poststim.png', subjid)), 'Resolution', 200);
    close(h);
end

% ================================================================
% FIGURE 3: Phase heatmap + GA rose plot
% polarhistogram needs a PolarAxes: the tile's regular axes is replaced
% by a polaraxes in the same position.
% ================================================================
if hasPhase
    phi_mat = featChan.slow_phase;

    % GA phase = circular mean across channels per trial
    phi_ga = angle(mean(exp(1i * phi_mat), 1, 'omitnan'));
    phi_ga_valid = phi_ga(~isnan(phi_ga));

    panelH = max(400, 14*nChan + 80);
    h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [50 50 900 panelH]);
    tl = tiledlayout(h, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('sub-%03d Slow-alpha Phase at t = 0', subjid), ...
        'Interpreter', 'tex', 'FontSize', 11, 'FontWeight', 'bold');

    ax1 = nexttile(tl, 1);
    imagesc(ax1, phi_mat);
    axis(ax1, 'tight');
    set(ax1, 'YDir', 'normal');
    set(ax1, 'YTick', yt, 'YTickLabel', chanLabels(yt), 'FontSize', 9);
    xlabel(ax1, 'Trial', 'FontSize', 8);
    ylabel(ax1, 'Channel', 'FontSize', 8);
    caxis(ax1, [-pi pi]);
    colormap(ax1, hsv(256));
    cb = colorbar(ax1, 'Location','eastoutside');
    cb.Ticks = [-pi -pi/2 0 pi/2 pi];
    cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};

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

    exportgraphics(h, fullfile(outDir, sprintf('sub-%03d_heatmap_phase.png', subjid)), 'Resolution', 200);
    close(h);
end

end

% ================================================================
% Local: one channel x trial heatmap
% ================================================================
function draw_heatmap(ax, X, ttl, signed, chanLabels, yt)
imagesc(ax, X);
axis(ax, 'tight');
set(ax, 'YDir', 'normal');
set(ax, 'YTick', yt, 'YTickLabel', chanLabels(yt), 'FontSize', 8);
title(ax, ttl, 'Interpreter', 'tex', 'FontSize', 9);
xlabel(ax, 'Trial', 'FontSize', 8);
ylabel(ax, 'Channel', 'FontSize', 8);
if signed
    clim_abs = max(abs(X(isfinite(X))));
    if ~isempty(clim_abs) && clim_abs > 0, caxis(ax, [-clim_abs clim_abs]); end
    colormap(ax, spec_diverging_cmap());
else
    colormap(ax, parula);
end
colorbar(ax, 'Location', 'eastoutside');
end

% ================================================================
% Local: keep only item rows whose field exists in featChan
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