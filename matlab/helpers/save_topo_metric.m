function save_topo_metric(LOGS, subjid, EEG, v, label)
h = figure('Visible', 'off');
topoplot(v, EEG.chanlocs, 'electrodes', 'on');
title(sprintf('sub-%03d %s topoplot', subjid, label), 'Interpreter', 'none');
saveas(h, fullfile(LOGS, sprintf('sub-%03d_initrej_topo_%s.png', subjid, lower(label))));
close(h);
end