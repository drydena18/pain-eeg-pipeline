function spec_plot_interaction_summary(outPath, featChan, featGA, chanLabels, subjid, logf)
%% SPEC_PLOT_INTERACTION_SUMMARY
% Per-subject overview of pre-stim interaction metrics
% V 1.1.0
%
% V 1.1.0 changes vs V 1.0.0:
%   - Fixed "Parent input must be a polaraxes object" error: nexttile returns
%     a regular Axes; polarhistogram requires a PolarAxes.  The fix captures
%     the tile's Position, deletes the regular axes, then creates a polaraxes
%     in the same screen location.
%   - Fixed missing semicolon on hasFlag assignment (was printing to console).
%   - Fixed title typo: "RL_pre" -> "LR_pre" in panel [1,1].
%
% Produces a single multi-panel figure covering all new alpha interaction metrics.
%
% Layout (3 rows x 3 columns):
%   Row 1 – GA trial sequence (line plots):
%           [1,1] BI_pre + LR_pre (dual y-axis)
%           [1,2] CoG_pre with 10 Hz boundary marked
%           [1,3] ΔERD with zero reference
%   Row 2 – Channel x trial heatmaps (diverging colourmap):
%           [2,1] BI_pre [chan x trial]
%           [2,2] ΔERD   [chan x trial]
%           [2,3] CoG_pre [chan x trial]
%   Row 3 – Phase / p5_flag / BI_pre autocorrelation:
%           [3,1] Slow-alpha phase polar histogram (rose)
%           [3,2] p5_flag heatmap [chan x trial]
%           [3,3] BI_pre GA autocorrelation (lag 1-10)
%
% Rows 3 panels are each conditional – drawn only when the relevant data exists.
% Missing panels display a "not available" placeholder.
%
% Inputs:
%   - outPath       : full path for the saved PNG
%   - featChan      : struct of [nChar x nTr] fields from spectral_core
%   - featGA        : struct of [nTr x 1] fields (GA across channels)
%   - chanLabels    : cell/string array of channel labels
%   - subjid        : integer subject ID
%   - logf          : (optional) MATLAB log file handle

if nargin < 6, logf = []; end
if isstring(chanLabels), chanLabels = cellstr(chanLabels); end
 
nTr = numel(featGA.bi_pre);
nC_panels = 3;
nR_panels = 3;
 
nChan = numel(chanLabels);
step = 1;
if nChan > 40, step = 2; end
if nChan > 80, step = 4; end
yt = 1:step:nChan;
 
figW = 1800;
figH = max(1400, 300*nR_panels + 20*nChan);
 
h = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 figW figH]);
tl = tiledlayout(h, nR_panels, nC_panels, 'TileSpacing', 'compact', 'Padding', 'compact');
 
title(tl, sprintf('sub-%03d | Pre-Stimulus Alpha Interaction Metrics', subjid), ...
    'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
 
trials = (1:nTr)';
 
% ================================================================
% ROW 1 – GA trial sequence
% ================================================================
% [1,1] BI_pre + LR_pre on dual axes
ax = nexttile(tl, 1);
yyaxis(ax, 'left');
plot(ax, trials, featGA.bi_pre, '-o', 'MarkerSize', 3, 'Color', [0.13 0.47 0.71]);
yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
ylabel(ax, 'BI_{pre} [-1, +1]', 'Interpreter', 'tex');
ylim(ax, [-1.05 1.05]);
 
yyaxis(ax, 'right');
plot(ax, trials, featGA.lr_pre, '-s', 'MarkerSize', 3, 'Color', [0.84 0.15 0.16]);
yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
ylabel(ax, 'LR_{pre} (log ratio)', 'Interpreter', 'tex');
 
title(ax, 'BI_{pre} & LR_{pre} over trials', 'Interpreter', 'tex');   % V1.1.0: was RL_pre
xlabel(ax, 'Trial');
legend(ax, {'BI_{pre}', 'LR_{pre}'}, 'Interpreter', 'tex', 'Location', 'best');
grid(ax, 'on');
 
% [1,2] CoG_pre – weighted spectral centroid within alpha band
ax = nexttile(tl, 2);
plot(ax, trials, featGA.cog_pre, '-o', 'MarkerSize', 3, 'Color', [0.17 0.63 0.17]);
yline(ax, 10, '--k', '10 Hz', 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
ylabel(ax, 'CoG_{pre} (Hz)', 'Interpreter', 'tex');
xlabel(ax, 'Trial');
title(ax, 'CoG_{pre}: Pre-stim Center of Gravity', 'Interpreter', 'tex');
grid(ax, 'on');
 
% [1,3] ΔERD – fast vs slow alpha desynchronisation asymmetry
ax = nexttile(tl, 3);
plot(ax, trials, featGA.delta_erd, '-o', 'MarkerSize', 3, 'Color', [0.58 0.40 0.74]);
yline(ax, 0, '--k', 'Equal', 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
ylabel(ax, '\DeltaERD (ERD_{slow} - ERD_{fast})', 'Interpreter', 'tex');
xlabel(ax, 'Trial');
title(ax, '\DeltaERD: Sub-band ERD asymmetry', 'Interpreter', 'tex');
grid(ax, 'on');
 
% ================================================================
% ROW 2 – Channel x trial heatmaps
% ================================================================
heatmap_items = {
    4, 'bi_pre',    'BI_{pre} [chan \times trial]';
    5, 'delta_erd', '\DeltaERD [chan \times trial]';
    6, 'cog_pre',   'CoG_{pre} [chan \times trial]';
};
 
for r = 1:size(heatmap_items, 1)
    tile_idx = heatmap_items{r, 1};
    fn  = heatmap_items{r, 2};
    ttl = heatmap_items{r, 3};
 
    ax = nexttile(tl, tile_idx);
 
    if isfield(featChan, fn) && ~isempty(featChan.(fn))
        X = featChan.(fn);
        imagesc(ax, X);
        axis(ax, 'tight');
        set(ax, 'YDir', 'normal');
        set(ax, 'YTick', yt, 'YTickLabel', chanLabels(yt), 'FontSize', 9);
        xlabel(ax, 'Trial');
        ylabel(ax, 'Channel');
        colorbar(ax, 'eastoutside');
        colormap(ax, spec_diverging_cmap());
        clim_abs = max(abs(X(:)), [], 'omitnan');
        if clim_abs > 0
            clim(ax, [-clim_abs, clim_abs]);
        end
    else
        text(ax, 0.5, 0.5, 'not available', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
        axis(ax, 'off');
    end
 
    title(ax, ttl, 'Interpreter', 'tex');
end
 
% ================================================================
% ROW 3 – Phase rose / p5_flag / BI_pre autocorrelation
% ================================================================
 
% ---------------------------------------------------------------
% [3,1] Circular histogram of GA phase at stimulus onset
%
% FIX (V1.1.0): nexttile returns a regular Axes object, but
% polarhistogram requires a PolarAxes as its parent.  We capture
% the tile Position, delete the regular axes, then create a
% polaraxes in the same screen rectangle.
% ---------------------------------------------------------------
ax_placeholder = nexttile(tl, 7);
hasPhase = isfield(featGA, 'phase_slow_rad') && ...
           ~all(isnan(featGA.phase_slow_rad));
 
if hasPhase
    phi = featGA.phase_slow_rad(~isnan(featGA.phase_slow_rad));
 
    % Replace the regular axes with a polaraxes in the same tile position
    pos = ax_placeholder.Position;
    delete(ax_placeholder);
    ax = polaraxes(h, 'Position', pos);
 
    polarhistogram(ax, phi, 24, 'Normalization', 'probability', ...
        'FaceColor', [0.13 0.47 0.71], 'FaceAlpha', 0.6);
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir          = 'clockwise';
    title(ax, 'GA Slow-alpha Phase at t = 0', 'Interpreter', 'none');
else
    ax = ax_placeholder;   % keep the regular axes for the placeholder
    text(ax, 0.5, 0.5, 'phase not available', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
    axis(ax, 'off');
    title(ax, 'GA Slow-alpha Phase at t = 0', 'Interpreter', 'none');
end
 
% ---------------------------------------------------------------
% [3,2] p5_flag heatmap – binary map of unstable-power trials
% ---------------------------------------------------------------
ax = nexttile(tl, 8);
hasFlag = isfield(featChan, 'p5_flag') && ~isempty(featChan.p5_flag);   % V1.1.0: added semicolon
 
if hasFlag
    F = featChan.p5_flag;
    imagesc(ax, F);
    axis(ax, 'tight');
    set(ax, 'YDir', 'normal');
    set(ax, 'YTick', yt, 'YTickLabel', chanLabels(yt), 'FontSize', 9);
    colormap(ax, [1 1 1; 0.84 0.15 0.16]);   % white = ok, red = flagged
    clim(ax, [0 1]);
    cb = colorbar(ax, 'eastoutside');
    cb.Ticks = [0 1]; cb.TickLabels = {'ok', 'flagged'};
    xlabel(ax, 'Trial'); ylabel(ax, 'Channel');
    pct = 100 * mean(F(:));
    title(ax, sprintf('p5 Flag (%.1f%% cells flagged)', pct), 'Interpreter', 'none');
else
    text(ax, 0.5, 0.5, 'not available', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
    axis(ax, 'off');
    title(ax, 'p5 Flag', 'Interpreter', 'none');
end
 
% ---------------------------------------------------------------
% [3,3] BI_pre autocorrelation (lags 1-10)
% Slow decay -> rigid state; rapid decay -> flexible state
% ---------------------------------------------------------------
ax = nexttile(tl, 9);
b    = featGA.bi_pre;
b_ok = b(~isnan(b));
maxLag = min(10, floor(numel(b_ok) / 3));
 
if numel(b_ok) >= 6
    [ac, lags] = xcorr(b_ok - mean(b_ok), maxLag, 'coeff');
    posLags = lags >= 0;
    bar(ax, lags(posLags), ac(posLags), 'FaceColor', [0.13 0.47 0.71]);
    yline(ax, 0, 'k');
    ci = 1.96 / sqrt(numel(b_ok));   % 95% CI under H0 (white noise)
    yline(ax, ci,  '--r', '95% CI', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
    yline(ax, -ci, '--r');
    xlabel(ax, 'Lag (trials)');
    ylabel(ax, 'Autocorrelation');
    title(ax, 'BI_{pre} GA Autocorrelation', 'Interpreter', 'tex');
    xlim(ax, [-0.5, maxLag + 0.5]);
    grid(ax, 'on');
else
    text(ax, 0.5, 0.5, 'insufficient trials', 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
    axis(ax, 'off');
    title(ax, 'BI_{pre} GA Autocorrelation', 'Interpreter', 'tex');
end
 
exportgraphics(h, outPath, 'Resolution', 200);
close(h);
 
spec_logmsg(logf, '[INTERACT_FIG] Saved -> %s', outPath);
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