function [phaseMat, rclTable] = spec_compute_phase_metric(EEG, alpha, phCfg, logf, ratings)
% SPEC_COMPUTE_PHASE_METRIC Slow-alpha Hilbert phase at stimulus onset
% V 1.1.0
%
% Computes the instantaneous phase of the slow-alpha [8, 10] Hz oscillation
% at t = 0 (stimulus onset) for every channel and trial. Phase is extracted
% from the analytic signal obtained by applying the Hilbert transform to a 
% zero-phase FIR-filtered copy of the epoched data.
%
% V 1.1.0 changes
%   - Accepts an explicit 'ratings' argument (5th parameter). If supplied
%     these are used for the circular-linear (r_cl) computation instead of
%     attempting to extract ratings from EEG.epoch.
%   - r_cl computation delegated to spec_compute_rcl_from_phase so that
%     the same logic can be used by the 09_hilbert loading path in
%     spectral_core without re-running the Hilbert transform.
%   - extract_epoch_ratings retained as a last-resort fallback when
%     ratings is [] and EEG.epoch.pain_rating is populated
%
% Inputs:
%   EEG         : EEGLAB epoched struct (preferably 08_base)
%   alpha       : struct with slow_hz [1 x 2]
%   phCfg       : cfg.spectral.phase: n_permutations (default 5000)
%   logf        : log file handle for spec_logmsg
%   ratings     : [nTrials x 1] pain_ratings loaded externally.
%
% Outputs:
%   phaseMat    : [nChan x nTrials] instantaneous phase at t=0 (radians, -pi to pi)
%   rclTable    : table(chan_idx, chan_label, r_cl, p_perm)

if nargin < 4, logf = 1; end
if nargin < 5, ratings = []; end

fs = EEG.srate;
nChan = EEG.nbchan;
nTr = EEG.trials;
slowHz = alpha.slow_hz;

% Stimulus-onset sample: the time point closest to 0 ms
timesSec = double(EEG.times(:))' / 1000;
[~, t0idx] = min(abs(timesSec));

spec_logmsg(logf, '[PHASE] Stimulus onset sample: id = %d (t = %.3f ms)', ...
    t0idx, EEG.times(t0idx));

% --------------------------------------------------
% Zero-phase FIR bandpass to slow alpha [8, 10] Hz
% pop_eegfiltnew pads and applies forward-backward FIR; it
% automatically enforces a sufficiently long filter order.
% Minimum order criterion: >= 3 * (fs / bandwidth)
% --------------------------------------------------
bw = slowHz(2) - slowHz(1);
minOrder = 3 * round(fs / bw);

spec_logmsg(logf, '[PHASE] Bandpass [%.0f, %.0f] Hz     FIR order >= %d     (fs = %.0f Hz)', ...
    slowHz(1), slowHz(2), minOrder, fs);

EEGfilt = pop_eegfiltnew(EEG, slowHz(1), slowHz(2));
EEGfilt = eeg_checkset(EEGfilt);

% --------------------------------------------------
% Hilbert transform per channel x trial
% Hilbert() operates on column vectors; pass [nTime x 1]
% --------------------------------------------------
phaseMat = nan(nChan, nTr);

for t = 1:nTr
    X = double(EEGfilt.data(:, :, t)); % [nChan x nTime]
    for ch = 1:nChan
        z = hilbert(X(ch, :)'); % [nTime x 1] complex analytic signal
        phaseMat(ch, t) = angle(z(t0idx));
    end
end

spec_logmsg(logf, '[PHASE] Phase extracted: %d chans x %d trials', nChan, nTr);

% --------------------------------------------------
% Resolve ratings: use explicitly passed ratings first
% fall back to EEG.epoch.pain_rating is ratings is empty
% --------------------------------------------------
if isempty(ratings)
    ratings = extract_epoch_ratings(EEG, nTr, logf);
end

% --------------------------------------------------
% Circular-linear correlation (delegated to shared helper)
% --------------------------------------------------
chanLabels = spec_get_chanlabels(EEG);
rclTable = spec_compute_rcl_from_phase(phaseMat, ratings, chanLabels, phCfg, logf);
end

% ==========================================================
% Local fallback: pull pain_rating out of EEG.epoch fields
% Used only when no external ratings vector is supplied
% ==========================================================
function ratings = extract_epoch_ratings(EEG, nTr, logf)
ratings = [];
if ~isfield(EEG, 'epoch') || isempty(EEG.epoch, 'pain_rating')
    return;
end
vals = nan(nTr, 1);
for t = 1:nTr
    ep = EEG.epoch(t);
    if isfield(ep, 'pain_rating')
        v = ep.pain_rating;
        if iscell(v), v = v{1}; end
        vals(t) = double(v);
    end
end
    if all(isnan(vals))
        spec_logmsg(logf, '[PHASE] EEG.epoch.pain_rating all NaN; no ratings available.');
        return;
    end
    if sum(~isnan(vals)) ~= nTr
        spec_logmsg(logf, '[PHASE][WARN] Only %d / %d trials have pain_rating in EEG.epoch', ...
            sum(~isnan(vals)), nTr);
    end
ratings = vals;
end