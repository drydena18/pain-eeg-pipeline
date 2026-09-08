function save_ic_topo_grid(QC, subjid, EEG, nICs)
% SAVE_IC_TOPO_GRID Save one figure with topomaps for the first N ICs.
%   Complements save_ic_qc_packets.m, which only generates full 4-panel
%   figures for ICLabel-flagged ICs. This gives a quick-look grid across
%   the first N components regardless of whether they were flagged, so
%   nothing outside ICLabel's suggestions goes unseen.
%
% Config: cfg.preproc.ica.grid.enabled / n_ics (default 20)
%
% No-op if EEG has no ICA decomposition or no channel locations.

if nargin < 4 || isempty(nICs)
    nICs = 20;
end

if ~isfield(QC, 'icawinv') || isempty(EEG.icawinv) || ~has_chanlocs(EEG)
    return;
end

outDir = fullfile(QC, sprintf('sub-%03d_icqc', subjid));
ensure_dir(outDir);

nTotal = size(EEG.icawinv, 2);
nShow = min(nICs, nTotal);

nCols = ceil(sqrt(nShow));
nRows = ceil(nShow / nCols);

C = [];
if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification') && isfield(EEG.etc.ica_classification, 'ICLabel') && isfield(EEG.etc.ica_classification.ICLabel, 'classifications')
    C = EEG.etc.ic_classification.ICLabel.classifications;
end
classNames = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line', 'ChanNoise', 'Other'};

h = figure('Visible', 'off', 'Position', [100 100 220 * nCols 200 * nRows]);

for ic = 1:nShow
    subplot(nRows, nCols, ic);
    topoplot(EEG.icawinv(:, ic), EEG.chanlocs, 'electrodes', 'off');

    ttl = sprintf('IC%d', ic);
    if ~isempty(C) && size(C, 1) >= ic
        [~, mi] = max(C(ic, :));
        ttl = sprintf('IC%d %s', ic, classNames{mi});
    end
    title(ttl, 'FontSize', 8, 'Interpreter', 'none');
end

sgtitle(sprintf('sub-%03d first %d ICs', subjid, nShow), 'Interpreter', 'none');

outPath = fullfile(outDir, sprintf('sub-%03d_ic_grid_first%d.png', subjid, nShow));
saveas(h, outPath);
close(h);

fprintf(1, '[ICQC] IC topomap grid witten to: %s\n', outPath);
end