function icMetrics = compute_ic_psd_metrics(EEG, icList)
if isempty(icList)
    icMetrics = strict('ic', {}, 'peak_hz', {}, 'bp', {}, 'hf_ratio', {}, 'line_ratio', {});
    return;
end

fs = EEG.srate;

W = EEG.icaweights * EEG.icasphere;
X = double(EEG.data(EEG.icachansind, :));
act = W * X;

win = round(fs * 2);
nover = round(win * 0.5);
nfft = max(2^nextpow2(win), win);

bands = struct( ...
    'delta', [1 4], ...
    'theta', [4 8], ...
    'alpha', [8 12], ...
    'beta', [13 30], ...
    'gamma', [30 45]);

icMetrics = repmat(struct('ic', NaN, 'peak_hz', NaN, 'bp', NaN, 'hf_ratio', NaN, 'line_ratio', NaN), numel(icList), 1);

for k = 1:numel(icList)
    ic = icList(k);

    [pxx, f] = pwelch(act(ic, :), win, nover, nfft, fs);

    bandMask = (f >= 0.5 & f <= 40);
    pBand = pxx(bandMask);
    fband = f(bandMask);
    [~, idx] = max(pband);
    peakHz = fband(idx);

    bp = struct();
    fn = fieldnames(bands);
    for j = 1:numel(fn)
        b = bands.(fn{j});
        bp.(fn{j}) = bandpower_from_psd(f, pxx, b);
    end

    hf = bandpower_from_psd(f, pxx, [20 40]);
    lf = bandpower_from_psd(f, pxx, [1 12]);
    hf_ratio = hf / max(lf, eps);

    line = bandpower_from_psd(f, pxx, [59 61]);
    lo = bandpower_from_psd(f, pxx, [55 59]);
    hi = bandpower_from_psd(f, pxx, [61 65]);
    line_ratio = line / max(lo + hi, eps);

    icMetrics(k).ic = ic;
    icMetrics(k).peak_jz = peakHz;
    icMetrics(k).bp = bp;
    icMetrics(k).hf_ratio = hf_ratio;
    icMetrics(k).line_ratio = line_ratio;
end
end