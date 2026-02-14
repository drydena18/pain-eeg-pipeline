function save_channel_psd_overview(LOGS, subjid, EEG)
fs = EEG.srate;
data = double(EEG.data);
nChan = size(data, 1);

win = round(fs * 2);
nover = round(win * 0.5);
nfft = max(2^nextpow2(win), win);

Pxx = zeros(nfft/2+1, nChan);
for ch = 1:nChan
    [pxx, f] = pwelch(data(ch, :), win, nover, nfft, fs);
    Pxx(:, ch) = 10*log10(pxx);
end

med = median(Pxx, 2);
q1 = prctile(Pxx, 25, 2);
q3 = prctile(Pxx, 75, 2);

h = figure('Visible', 'off');
plot(f, med); hold on;
plot(f, q1, ':'); plot(f, q3, ':');
xlim([0 80]);
xlabel('Hz'); ylabel('Power (dB)');
title(sprintf('sub-%03d channel PSD (median + IQR)', subjid), 'Interpreter', 'none');
legend({'Median', '25%', '75%'}, 'Location', 'northeast');
saveas(h, fullfile(LOGS, sprintf('sub-%03d_initrej_psd_overview.png', subjid)));
close(h);
end