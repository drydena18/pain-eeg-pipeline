function save_channel_psd_badchans(LOGS, subjid, EEG, badChans)
fs = EEG.srate;
data = double(EEG.data);

win = round(fs * 2);
nover = round(win * 0.5);
nfft = max(2^nextpow2(win), win);

h = figure('Visible', 'off'); hold on;
for ch = badChans(:)'
    [pxx, f] = pwelch(data(ch, :), win, nover, nfft, fs);
    plot(f, 10*log10(pxx));
end
xlim([0 80]);
xlabel('Hz'); ylabel('Power (dB)');
title(sprintf('sub-%03d PSD overlay (suggested channels %s)', subjid, vec2str(badChans)), 'Interpreter', 'none');
saveas(h, fullfile(LOGS, sprintf('sub-%03d_initrej_psd_badchans.png', subjid)));
close(h);
end