function save_hist(LOGS, subjid, v, name)
h = figure('Visible', 'off');
histogram(v);
title(sprintf('sub-%03d %s histogram', subjid, name), 'Interpreter', 'none');
xlabel(name); ylabel('Count');
saveas(h, fullfile(LOGS, sprintf('sub-%03d_initrej_hist_%s.png', subjid, name)));
close(h);
end