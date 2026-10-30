function [EEG, qc] = st_hilbert(EEG, cfg, ~)
%ST_HILBERT Stage 09 (preproc_core slow-alpha Hilbert block) ->
%   09_hilbert/sub-XXX_hilbert_phase.mat. EEG is not modified. Unlike
%   preproc_core, an existing .mat is overwritten: Mode 2b only runs when
%   the IC decisions or post-ICA settings changed.
qc = struct();
h = struct('enabled', false);
if isfield(cfg.exp.preproc, 'hilbert'), h = cfg.exp.preproc.hilbert; end
if ~(isfield(h, 'enabled') && logical(h.enabled)), qc.disabled = true; return; end

rt = cfg.rt; logf = rt.logf; subjid = rt.subjid; ecfg = cfg.exp;
matPath = fullfile(rt.ST.HILBERT, sprintf('sub-%03d_hilbert_phase.mat', subjid));

if EEG.trials <= 1
    logmsg(logf, '[HILBERT][WARN] EEG is not epoched (trials = %d); stage 09 skipped.', EEG.trials);
    qc.skipped = true;
    return;
end

slowHz = [8 10];
if isfield(ecfg, 'features') && isfield(ecfg.features, 'bands') && isfield(ecfg.features.bands, 'slow_hz')
    slowHz = double(ecfg.features.bands.slow_hz(:)');
elseif isfield(ecfg, 'spectral') && isfield(ecfg.spectral, 'alpha') && isfield(ecfg.spectral.alpha, 'slow_hz')
    slowHz = ecfg.spectral.alpha.slow_hz;
elseif isfield(h, 'slow_hz') && ~isempty(h.slow_hz)
    slowHz = h.slow_hz;
end

fs = EEG.srate; nChan = EEG.nbchan; nTr = EEG.trials; nTime = EEG.pnts;
timesSec = double(EEG.times(:))' / 1000;
[~, t0idx] = min(abs(timesSec));
logmsg(logf, '[HILBERT] Band = [%.0f %.0f] Hz nChan = %d nTrials = %d stimOnsetIdx = %d', ...
    slowHz(1), slowHz(2), nChan, nTr, t0idx);

EEGfilt = pop_eegfiltnew(EEG, slowHz(1), slowHz(2));
EEGfilt = eeg_checkset(EEGfilt);

phaseSlow = zeros(nChan, nTime, nTr, 'single');
for t = 1:nTr
    X = double(EEGfilt.data(:, :, t));
    for ch = 1:nChan
        phaseSlow(ch, :, t) = single(angle(hilbert(X(ch, :)')));
    end
end

phase_slow   = phaseSlow; %#ok<NASGU>
t_ms         = EEG.times; %#ok<NASGU>
stimOnsetIdx = t0idx;     %#ok<NASGU>
chan_labels  = get_chan_labels(EEG, sprintf('sub-%03d stage09', subjid)); %#ok<NASGU>
slow_hz      = slowHz;    %#ok<NASGU>
subjid_save  = subjid;    %#ok<NASGU>
save(matPath, 'phase_slow', 't_ms', 'stimOnsetIdx', 'chan_labels', 'slow_hz', 'subjid_save', 'fs', '-v7.3');
logmsg(logf, '[HILBERT] Saved -> %s', matPath);

qc.mat = matPath;
qc.slow_hz = slowHz;
end