function [phaseMat, rclTable] = spec_compute_phase_metric(EEG, alpha, phCfg, logf)
% SPEC_COMPUTE_PHASE_METRIC  Slow-alpha Hilbert phase at stimulus onset
% V 1.0.0
%
% Computes the instantaneous phase of the slow-alpha [8,10] Hz oscillation
% at t = 0 (stimulus onset) for every channel and trial.  Phase is extracted
% from the analytic signal obtained by applying the Hilbert transform to a
% zero-phase FIR-filtered copy of the epoched data.
%
% Optionally computes the circular-linear correlation (Harrison-Kanji r_cl)
% between phase and pain rating per channel, with a permutation p-value.
% This requires pain_rating to be embedded in EEG.epoch (field 'pain_rating').
%
% If pain ratings are unavailable, phaseMat is still returned and rclTable
% is an empty table.  The modelling in R uses sin/cos of phase directly.
%
% Inputs:
%   EEG    : EEGLAB epoched struct (preferably 08_base)
%   alpha  : struct with slow_hz [1x2]  e.g. [8, 10]
%   phCfg  : cfg.spectral.phase:
%               n_permutations  (default 5000)
%   logf   : (optional) MATLAB log file handle
%
% Outputs:
%   phaseMat : [nChan x nTrials]  instantaneous phase at t=0 (radians, -pi to pi)
%   rclTable : table with chan_idx, chan_label, r_cl, p_perm  (empty if no ratings)

if nargin < 4, logf = 1; end

nPerm    = 5000;
if isfield(phCfg, 'n_permutations') && ~isempty(phCfg.n_permutations)
    nPerm = phCfg.n_permutations;
end

fs      = EEG.srate;
nChan   = EEG.nbchan;
nTr     = EEG.trials;
slowHz  = alpha.slow_hz;   % [8 10]

% Stimulus-onset sample: the time point closest to 0 ms
timesSec  = double(EEG.times(:))' / 1000;
[~, t0idx] = min(abs(timesSec));

spec_logmsg(logf, '[PHASE] Stimulus onset sample: idx=%d (t=%.3f ms)', ...
    t0idx, EEG.times(t0idx));

% ---------------------------------------------------------------
% Zero-phase FIR bandpass to slow alpha [8, 10] Hz
% EEGLAB's pop_eegfiltnew pads and applies forward-backward FIR.
% Minimum order criterion from the markdown: >= 3 * (fs / df)
% where df = bandwidth = slow_hz(2) - slow_hz(1).
% ---------------------------------------------------------------
bw        = slowHz(2) - slowHz(1);
minOrder  = 3 * round(fs / bw);

spec_logmsg(logf, '[PHASE] Bandpass [%.0f %.0f] Hz  FIR order >= %d  (fs=%.0f Hz)', ...
    slowHz(1), slowHz(2), minOrder, fs);

% pop_eegfiltnew(EEG, locutoff, hicutoff) applies a zero-phase FIR.
% The minimum filter order is controlled internally by EEGLAB to be at
% least (3 * fs / locutoff) cycles, which is conservative enough.
EEGfilt = pop_eegfiltnew(EEG, slowHz(1), slowHz(2));
EEGfilt = eeg_checkset(EEGfilt);

% ---------------------------------------------------------------
% Hilbert transform per channel x trial
% MATLAB hilbert() operates on column vectors; we pass [nTime x 1].
% ---------------------------------------------------------------
phaseMat = nan(nChan, nTr);

for t = 1:nTr
    X = double(EEGfilt.data(:, :, t));  % [nChan x nTime]
    for ch = 1:nChan
        z                = hilbert(X(ch, :)');    % analytic signal [nTime x 1]
        phaseMat(ch, t)  = angle(z(t0idx));       % instantaneous phase at t=0
    end
end

spec_logmsg(logf, '[PHASE] Phase extracted: %d chans x %d trials.', nChan, nTr);

% ---------------------------------------------------------------
% Circular-linear correlation (optional; requires ratings in EEG.epoch)
% ---------------------------------------------------------------
rclTable = table();

ratings = extract_epoch_ratings(EEG, nTr, logf);
if isempty(ratings)
    spec_logmsg(logf, '[PHASE] No pain ratings in EEG.epoch; skipping r_cl computation.');
    return;
end

chanLabels = spec_get_chanlabels(EEG);
r_cl_vec   = nan(nChan, 1);
p_perm_vec = nan(nChan, 1);

for ch = 1:nChan
    phi = phaseMat(ch, :)';    % [nTr x 1]
    ok  = ~isnan(phi) & ~isnan(ratings);

    if sum(ok) < 10
        spec_logmsg(logf, '[PHASE][WARN] Chan %d: only %d valid trials; skipping r_cl.', ch, sum(ok));
        continue;
    end

    r_cl_vec(ch)   = circ_lin_corr(phi(ok), ratings(ok));
    p_perm_vec(ch) = permutation_p(phi(ok), ratings(ok), nPerm);
end

rclTable = table( ...
    (1:nChan)', string(chanLabels(:)), r_cl_vec, p_perm_vec, ...
    'VariableNames', {'chan_idx', 'chan_label', 'r_cl', 'p_perm'} ...
);

spec_logmsg(logf, '[PHASE] r_cl range [%.3f, %.3f] (median p_perm = %.3f)', ...
    min(r_cl_vec, [], 'omitnan'), ...
    max(r_cl_vec, [], 'omitnan'), ...
    median(p_perm_vec, 'omitnan'));
end

% ================================================================
%  Local: Harrison-Kanji circular-linear correlation coefficient
%  r_cl = sqrt( (r_cs^2 + r_cc^2 - 2*r_cs*r_cc*r_sc) / (1 - r_sc^2) )
% ================================================================
function r = circ_lin_corr(phi, y)
s    = sin(phi);
c    = cos(phi);
r_cs = corr(y, s, 'type', 'Pearson');
r_cc = corr(y, c, 'type', 'Pearson');
r_sc = corr(s,  c, 'type', 'Pearson');
denom = 1 - r_sc^2;
if abs(denom) < eps
    r = 0; return;
end
r2 = (r_cs^2 + r_cc^2 - 2*r_cs*r_cc*r_sc) / denom;
r  = sqrt(max(r2, 0));
end

% ================================================================
%  Local: permutation p-value (one-sided upper tail, H0: r_cl = 0)
%  Uses (count + 1) / (nPerm + 1) to avoid p = 0 (Phipson & Smyth 2010)
% ================================================================
function p = permutation_p(phi, y, nPerm)
obs   = circ_lin_corr(phi, y);
nT    = numel(phi);
count = 0;
for k = 1:nPerm
    phi_s = phi(randperm(nT));
    count = count + (circ_lin_corr(phi_s, y) >= obs);
end
p = (count + 1) / (nPerm + 1);
end

% ================================================================
%  Local: try to pull pain_rating out of EEG.epoch fields
% ================================================================
function ratings = extract_epoch_ratings(EEG, nTr, logf)
ratings = [];
if ~isfield(EEG, 'epoch') || isempty(EEG.epoch)
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
    spec_logmsg(logf, '[PHASE][WARN] Only %d / %d trials have pain_rating.', sum(~isnan(vals)), nTr);
end
ratings = vals;
end