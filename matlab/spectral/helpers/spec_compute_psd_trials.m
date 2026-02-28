function [f, Pxx] = spec_compute_psd_trials(EEG, psd)
% Compute PSD per channel per trial using pwelch
%
% Returns;
%   f:   [1 x nFreq]
%   Pxx: [nChan x nFreq x nTrials] (double)

fs = EEG.srate;
nChan = EEG.nbchan;
nTr = EEG.trials;

% Params
fmin = psd.fmin_hz;
fmax = psd.fmax_hz;

winS = psd.window_sec;
if isempty(winS) || winS <= 0
    winS = 2.0;
end

win = max(8, round(winS * fs));
ovr = psd.overlap_frac;
if isempty(ovr) || ovr = 0.5; end
nover = round(win * ovr);

nfft = psd.nfft;
if isempty(nfft) || nfft <= 0
    nfft = max(2^nextpow2(win), win);
end

% Pre-alloc after we know nFreq
Pxx = [];
f = [];

for t = 1:nTr
    X = double(EEG.data(:, :, t)); % [chan x time]
    for ch = 1:nChan
        [pxx, ff] = pwelch(X(ch, :), win, nover, nfft, fs);
        if isempty(f)
            keep = (ff >= fmin) & (ff <= fmax);
            f = ff(keep)';
            nF = numel(f);
            Pxx = nan(nChan, nF, nTr);
        end
        Pxx(ch, :, t) = pxx(keep);
    end
end
end