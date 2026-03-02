function spec_plot_trial_spectral_qc(outDir, EEG, cfg, subjid, logf)
% Saves one figure per trial: Pre-stim vs post-stim PSD (all channels)
% Uses a COMMON Welch window/nfft for pre+post so frequency grids match.

ts = cfg.spectral.trial_spectral;

preSec   = ts.pre_sec;
postSec  = ts.post_sec;
maxTrials = ts.max_trials;
fmin = ts.fmin_hz;
fmax = ts.fmax_hz;

psd = cfg.spectral.psd; % uses window_sec/overlap_frac/nfft as "desired"

legendMax = 16;
if isfield(ts, 'legend_max_channels'), legendMax = ts.legend_max_channels; end

if EEG.trials <= 1
    spec_logmsg(logf, '[TRIALSPEC][WARN] EEG not epoched. Skipping...');
    return;
end

if ~isfield(EEG, 'times') || isempty(EEG.times)
    error('spec_plot_trial_spectral_qc:MissingTimes', 'EEG.times is missing; cannot split pre/post.');
end

timesSec = double(EEG.times(:))' / 1000; % seconds
idxPre  = timesSec >= preSec(1)  & timesSec <= preSec(2);
idxPost = timesSec >= postSec(1) & timesSec <= postSec(2);

nPre  = nnz(idxPre);
nPost = nnz(idxPost);

if nPre < 8 || nPost < 8
    spec_logmsg(logf, '[TRIALSPEC][WARN] Not enough samples in pre/post windows. pre=%d post=%d', nPre, nPost);
    return;
end

fs = EEG.srate;

% ---- Desired Welch window from cfg ----
winS = psd.window_sec;
if isempty(winS) || winS <= 0, winS = 2.0; end
win_desired = max(8, round(winS * fs));

ovr = psd.overlap_frac;
if isempty(ovr) || ovr <= 0 || ovr >= 1, ovr = 0.5; end

nfft_cfg = psd.nfft;
if isempty(nfft_cfg) || nfft_cfg <= 0, nfft_cfg = 0; end

% ---- COMMON window based on the SHORTER segment (pre usually) ----
win = min(win_desired, min(nPre, nPost));
win = max(win, 8);

nover = round(win * ovr);
if nover >= win, nover = max(0, win - 1); end

if nfft_cfg > 0
    nfft = max(nfft_cfg, win);
else
    nfft = max(2^nextpow2(win), win);
end

spec_logmsg(logf, '[TRIALSPEC] Welch common params: win=%d samp (%.3fs), nover=%d, nfft=%d', ...
    win, win/fs, nover, nfft);

chanLabels = spec_get_chanlabels(EEG);
nChan = EEG.nbchan;
nTr = EEG.trials;

nDo = min(nTr, maxTrials);
if nDo < nTr
    spec_logmsg(logf, '[TRIALSPEC] Limiting trials: %d/%d', nDo, nTr);
end

for t = 1:nDo
    try
        Xpre  = double(EEG.data(:, idxPre,  t));  % [chan x nPre]
        Xpost = double(EEG.data(:, idxPost, t));  % [chan x nPost]

        % Prealloc based on first channel
        [pxx0, ff] = pwelch(Xpre(1,:), win, nover, nfft, fs);
        keep = (ff >= fmin) & (ff <= fmax);
        f = ff(keep)';

        pPre  = nan(nChan, numel(f));
        pPost = nan(nChan, numel(f));

        pPre(1,:) = pxx0(keep);

        % Remaining channels
        for ch = 2:nChan
            pxx = pwelch(Xpre(ch,:), win, nover, nfft, fs);
            pPre(ch,:) = pxx(keep);
        end

        for ch = 1:nChan
            pxx = pwelch(Xpost(ch,:), win, nover, nfft, fs);
            pPost(ch,:) = pxx(keep);
        end

        outPath = fullfile(outDir, sprintf('sub-%03d_trial-%03d_prepost_psd.png', subjid, t));
        spec_plot_trial_prepost_psd(outPath, f, pPre, f, pPost, chanLabels, subjid, t, legendMax);

    catch ME
        spec_logmsg(logf, '[TRIALSPEC][WARN] Trial %d failed: %s', t, ME.message);
    end
end
end