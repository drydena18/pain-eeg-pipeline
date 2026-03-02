function spec_plot_trial_spectral_qc(outDir, EEG, cfg, subjid, logf)
% Saves one figure per trial: Pre-stim vs post-stim PSD (all channels)
%
% Requires epoched EEG with EEG.times in ms.

% ---- Pull config ----
ts = cfg.spectral.trial_spectral;

preSec    = ts.pre_sec;
postSec   = ts.post_sec;
maxTrials = ts.max_trials;

fmin = ts.fmin_hz;
fmax = ts.fmax_hz;

% Use same pwelch params as cfg.spectral.psd
psd = cfg.spectral.psd;

legendMax = 16;
if isfield(ts, 'legend_max_channels') && ~isempty(ts.legend_max_channels)
    legendMax = ts.legend_max_channels;
end

% ---- Guardrails ----
if EEG.trials <= 1
    spec_logmsg(logf, '[TRIALSPEC][WARN] EEG not epoched. Skipping...');
    return;
end

if ~isfield(EEG, 'times') || isempty(EEG.times)
    error('spec_plot_trial_spectral_qc:MissingTimes', 'EEG.times is missing; cannot split pre/post.');
end

timesSec = double(EEG.times(:))' / 1000; % ms -> sec

idxPre  = timesSec >= preSec(1)  & timesSec <= preSec(2);
idxPost = timesSec >= postSec(1) & timesSec <= postSec(2);

if nnz(idxPre) < 8 || nnz(idxPost) < 8
    spec_logmsg(logf, '[TRIALSPEC][WARN] Not enough samples in pre/post windows. pre=%d post=%d', nnz(idxPre), nnz(idxPost));
    return;
end

chanLabels = spec_get_chanlabels(EEG);
nTr = EEG.trials;

nDo = min(nTr, maxTrials);
if nDo < nTr
    spec_logmsg(logf, '[TRIALSPEC] Limiting trials: %d/%d', nDo, nTr);
end

spec_logmsg(logf, '[TRIALSPEC] Writing %d pre/post PSD plots to: %s', nDo, outDir);

% ---- Main loop ----
for t = 1:nDo
    try
        [fPre,  pPre]  = spec_compute_psd_epoch_window(EEG, psd, t, idxPre,  fmin, fmax);
        [fPost, pPost] = spec_compute_psd_epoch_window(EEG, psd, t, idxPost, fmin, fmax);

        outPath = fullfile(outDir, sprintf('sub-%03d_trial-%03d_prepost_psd.png', subjid, t));
        spec_plot_trial_prepost_psd(outPath, fPre, pPre, fPost, pPost, chanLabels, subjid, t, legendMax);

    catch ME
        spec_logmsg(logf, '[TRIALSPEC][WARN] Trial %d failed: %s', t, ME.message);
    end
end

end