function chanPSD = compute_channel_psd_metrics(EEG)
fs = EEG.srate;
data = double(EEG.data);
nChan = size(data, 1);

win = round(fs * 2);
nover = round(win * 0.5);
nfft = max(2^nextpow2(win), win);

chanPSD = struct();
chanPSD.line_ratio  = zeros(nChan, 1);
chanPSD.hf_ratio    = zeros(nChan, 1);
chanPSD.drift_ratio = zeros(nChan, 1);
chanPSD.alpha_ratio = zeros(nChan, 1);

for ch = 1:nChan
    [pxx, f] = pwelch(data(ch, :), win, nover, nfft, fs);

    line = bp_psd(f, pxx, [59 61]);
    lo = bp_psd(f, pxx, [55 59]);
    hi = bp_psd(f, pxx, [61 65]);
    chanPSD.line_ratio(ch) = line / max(lo + hi, eps);

    hf = bp_psd(f, pxx, [20 40]);
    lf = bp_psd(f, pxx, [1 12]);
    chanPSD.hf_ratio(ch) = hf / max(lf, eps);

    drift = bp_psd(f, pxx, [1 2]);
    chanPSD.drift_ratio(ch) = drift / max(lf, eps);

    a = bp_psd(f, pxx, [8 12]);
    allp = bp_psd(f, pxx, [1 40]);
    chanPSD.alpha_ratio(ch) = a / max(allp, eps);
end
end