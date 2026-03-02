function [f, Pxx] = spec_compute_psd_epoch_window(EEG, psd, trialIdx, timeMask, fmin, fmax)
% Compute PSD for one trial and one time window mask (pre or post)
% Returns:
%   f   : [1 x nFreq]
%   Pxx : [nChan x nFreq]  (power)

fs = EEG.srate;

X = double(EEG.data(:, timeMask, trialIdx)); % [chan x nSamp]
nChan = size(X, 1);
nSamp = size(X, 2);

if nSamp < 8
    error('spec_compute_psd_epoch_window:TooShort', ...
        'Window too short: %d samples.', nSamp);
end

% --- Desired Welch params from cfg ---
winS = psd.window_sec;
if isempty(winS) || winS <= 0
    winS = 2.0;
end
win_desired = max(8, round(winS * fs));

ovr = psd.overlap_frac;
if isempty(ovr) || ovr <= 0 || ovr >= 1
    ovr = 0.5;
end

nfft_cfg = psd.nfft;
if isempty(nfft_cfg) || nfft_cfg <= 0
    nfft_cfg = 0; % use auto
end

% --- Adapt window to segment length ---
win = min(win_desired, nSamp);          % MUST be <= nSamp
if win < 8
    win = nSamp;                         % last resort
end

nover = round(win * ovr);
if nover >= win
    nover = max(0, win - 1);
end

% nfft must be >= win
if nfft_cfg > 0
    nfft = max(nfft_cfg, win);
else
    nfft = max(2^nextpow2(win), win);
end

% Compute for first channel to define f/keep
[pxx0, ff] = pwelch(X(1, :), win, nover, nfft, fs);
keep = (ff >= fmin) & (ff <= fmax);
f = ff(keep)';

nF = numel(f);
Pxx = nan(nChan, nF);

Pxx(1, :) = pxx0(keep);

for ch = 2:nChan
    pxx = pwelch(X(ch, :), win, nover, nfft, fs);
    Pxx(ch, :) = pxx(keep);
end
end