function spec_plot_debug_trials(outFirDir, f, Pxx, featChan, chanLabels, alpha, cfg,subjid, plotMode)
% Debug/exhaustive plots:
%   - Select subset of trials (debug) or all trials (exhaustive)
%   - Plot multi-channel PSD overlay wth limited legend

nTr = size(Pxx, 3);
qc = cfg.spectral.qc;

if plotMode == "exhaustive"
    trials = 1:nTr;
else
    trials = spec_select_trials_to_plot(featChan, qc, nTr);
end

legendMax = 20;
if isfield(qc, 'legend_max_channels'), legendMax = qc.legend_max_channels; end

for t = trials
    h = figure('Visible', 'off');
    hold on;
    for ch = 1:size(Pxx, 1)
        plot(f, 10*log10(Pxx(ch, :, t)));
    end
    xline(alpha.alpha_hz(1)); xline(alpha.alpha_hz(2));
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    title(sprintf('sub-%03d trial %d: PSD (all channels)', subjid, t), 'Interpreter', 'none');

    % Legend subset
    if numel(chanLabels) <= legendMax
        legend(cellstr(chanLabels), 'Interpreter', 'none', 'Location', 'eastoutside');
    else
        legend(cellstr(chanLabels(1:legendMax)), 'Interpreter', 'none', 'Location', 'eastoutside');
    end
    
    outPath = fullfile(outFigDir, sprintf('sub-%03d_trial%03d_psd_channels.png', subjid, t));
    exportgraphics(h, outPath, 'Resolution', 150);
    close(h);
end
end