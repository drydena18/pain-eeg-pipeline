function save_bar(LOGS, subjid, v, badChans, name)
h = figure('Visible', 'off');
bar(v);
title(sprintf('sub-%03d %s (suggested marked)', subjid, name), 'Interpreter', 'none');
xlabel('Channel index'); ylabel(name);
hold on;
for c = badChans(:)'
    xline(c, '--');
end
hold off;
saveas(h, fullfile(LOGS, sprintf('sub-%03d)initrej_bar_%s.png', subjid, name)));
close(h);
end