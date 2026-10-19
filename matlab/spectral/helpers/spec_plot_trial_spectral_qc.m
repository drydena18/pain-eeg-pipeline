function spec_plot_trial_spectral_qc(outDir, EEG, cfg, subjid, logf)
% SPEC_PLOT_TRIAL_SPECTRAL_QC
% V 1.1.0
%
% Saves one spectrogram figure per trial: Hz on Y, time (ms) on X,
% power (dB) on colour.  Channels are averaged in the linear power domain
% before plotting so the result is the grand-average (GA) spectrogram.
%
% V1.1.0: cfg.spectral.psd no longer exists (PSD settings moved to the
%         shared `features` block in V3), and reading it crashed the whole
%         QC step. The only value used from it was an optional nfft
%         override; it is now read from cfg.spectral.trial_spectral.nfft if
%         set, otherwise chosen automatically.
%
% Config keys read (all under cfg.spectral.trial_spectral):
%   max_trials             - cap on number of trials saved  (default 20)
%   fmin_hz / fmax_hz      - frequency display range        (default 1/40)
%   spectrogram_win_sec    - STFT window length in seconds  (default 0.5)
%   spectrogram_overlap    - overlap fraction 0-1           (default 0.90)
%   nfft                   - optional FFT length            (default auto)

ts = cfg.spectral.trial_spectral;

% ---- Basic guards -------------------------------------------------------
if EEG.trials <= 1
    spec_logmsg(logf, '[TRIALSPEC][WARN] EEG not epoched. Skipping.');
    return;
end

if ~isfield(EEG, 'times') || isempty(EEG.times)
    spec_logmsg(logf, '[TRIALSPEC][WARN] EEG.times missing. Skipping.');
    return;
end

% ---- Params -------------------------------------------------------------
maxTrials = defaultFieldVal(ts, 'max_trials', 20);
fmin      = defaultFieldVal(ts, 'fmin_hz',    1);
fmax      = defaultFieldVal(ts, 'fmax_hz',    40);

% Spectrogram-specific window: shorter than the PSD window to give useful
% time resolution.  0.5 s -> ~2 Hz freq resolution at typical EEG rates,
% which is sufficient to resolve the 8-10 and 10-12 Hz sub-bands.
winSec   = defaultFieldVal(ts, 'spectrogram_win_sec', 0.5);
ovrFrac  = defaultFieldVal(ts, 'spectrogram_overlap',  0.90);
nfft_cfg = defaultFieldVal(ts, 'nfft', 0);

fs   = EEG.srate;
nChan = EEG.nbchan;
nTr   = EEG.trials;
epochStartMs = double(EEG.times(1)); % ms relative to event (negative = pre-stim)

% Convert to samples
win   = max(8, round(winSec * fs));
nover = min(win - 1, round(win * ovrFrac));

if nfft_cfg > 0
    nfft = max(nfft_cfg, win);
else
    nfft = max(2^nextpow2(win), win);
end

spec_logmsg(logf, '[TRIALSPEC] Spectrogram: win=%d samp (%.3fs), noverlap=%d, nfft=%d, frange=[%.0f %.0f] Hz', ...
    win, win/fs, nover, nfft, fmin, fmax);

alpha = cfg.spectral.alpha;

nDo = min(nTr, maxTrials);
if nDo < nTr
    spec_logmsg(logf, '[TRIALSPEC] Capping at %d/%d trials.', nDo, nTr);
end

% ---- Per-trial loop -----------------------------------------------------
for t = 1:nDo
    try
        X = double(EEG.data(:, :, t)); % [nChan x nSamp]

        % Accumulate linear power across channels
        GA_S = [];
        f    = [];
        tMs  = [];

        for ch = 1:nChan
            [S, f_raw, t_raw] = spectrogram(X(ch, :), win, nover, nfft, fs);

            if ch == 1
                % Define frequency mask once
                fkeep = (f_raw >= fmin) & (f_raw <= fmax);
                f     = f_raw(fkeep);             % [nF x 1]
                % t_raw is seconds from first sample; shift to epoch ms
                tMs   = epochStartMs + t_raw(:)' * 1000; % [1 x nFrames]
                GA_S  = zeros(sum(fkeep), numel(t_raw));
            end

            GA_S = GA_S + abs(S(fkeep, :)).^2; % accumulate power
        end

        GA_S = GA_S / nChan; % channel average

        outPath = fullfile(outDir, ...
            sprintf('sub-%03d_trial-%03d_spectrogram.png', subjid, t));

        spec_plot_trial_spectrogram(outPath, GA_S, f, tMs, alpha, subjid, t);

    catch ME
        spec_logmsg(logf, '[TRIALSPEC][WARN] Trial %d failed: %s', t, ME.message);
    end
end

end % ---- end main -------------------------------------------------------

% ---- Local helper (avoids dependency on defaultField) -------------------
function v = defaultFieldVal(s, fn, def)
    if isfield(s, fn) && ~isempty(s.(fn))
        v = s.(fn);
    else
        v = def;
    end
end