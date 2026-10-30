function pp_plot_ic_grid(EEG, nShow, probs, classes, suggested, pngFile, subj)
%PP_PLOT_IC_GRID Topoplot grid of the first nShow ICs (runica orders them
%   by variance). Suggested rejections are titled in red.

nc = 10;
nr = ceil(nShow / nc);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 1700 175 * nr]);
locs = EEG.chanlocs(EEG.icachansind);

for k = 1:nShow
    subplot(nr, nc, k);
    topoplot(EEG.icawinv(:, k), locs, 'electrodes', 'off');
    [pk, jk] = max(probs(k, :));
    col = 'k';
    if ismember(k, suggested), col = 'r'; end
    title(sprintf('IC%d %s %.2f', k, classes{jk}, pk), 'Color', col, 'FontSize', 7);
end
sgtitle(sprintf('%s - %d of %d ICs (red = suggested reject)', subj, nShow, size(probs, 1)), ...
    'Interpreter', 'none');

exportgraphics(fig, pngFile, 'Resolution', 110);
close(fig);
end