function [fKeep, PxxKeep] = spec_compute_psd_epoch_window(EEG, psd, trialIdx, timeMask, fmin, fmax)
% Returns:
%   fKeep: [nFreq x 1] 
%   PxxKeep: [nChan x nFreq]

fs = EEG.srate;
X = double(EEG.data(:, timeMask, trialIdx));

% Welch params
winS = psd.window_sec;
if isempty(winS) || winS <= 0, winS = 2.0; end

win = max(8, round(winS * fs));
nover = round(win * psd.overlap_frac);

nfft = psd.nfft;
if isempty(nfft) || nfft <= 0
    nfft = max(2^nextpow2(win), win);
end

nChan = EEG.nbchan;

fKeep = [];
PxxKeep = [];

for ch = 1:nChan
    [pxx, ff] = pwelch(X(ch, :), win, nover, nfft, fs);

    if isempty(fKeep)
        keep = (ff >= fmin) & (ff <= fmax);
        fKeep = ff(keep)';
        PxxKeep = nan(nChan, numel(fKeep));
    end

    PxxKeep(ch, :) = pxx(keep);
end

end