function [featChan, featGA] = spec_compute_interaction_metrics(f, prePxx, postPxx, alpha, logf)
% SPEC_COMPUTE_INTERACTION_METRICS  Pre-stimulus alpha interaction metrics
% V 2.1.0
%
% V2.1.0: 
%   - Added explicit nargin / ndims / size-match check on prePxx and
%     postPxx BEFORE any band-power computation.
%
% V2.0.0 changes vs V1.1.0:
%   - REMOVED bi_pre, lr_pre, cog_pre, psi_cog. These are not computed
%     generically for every window (whole/pre/post) directly in spectral_core.m
%     via spec_compute_alpha_features_from_psd + spec_compute_psi_cog,
%     with delta_<name> from spec_compute_metric_deltas.m. This function now owns
%     ONLY the fractional (ERD-style) pre->post change, which is NOT produced
%     generically because it is not meaningful for bounded/ratio metrics
%   - ADDED erd_pow_alpha_total: fractional pre->post change for total
%     alpha power, completing the family alongside the pre-existing erd_slow / erd_fast.
%
% Computes the following per-channel x trial fields:
%   erd_slow   : Slow-alpha ERD  (post-pre)/pre  (negative = desynchronisation)
%   erd_fast   : Fast-alpha ERD  (post-pre)/pre
%   erd_pow_alpha_total : (post_alpha - pre_alpha) / pre_alpha
%   delta_erd  : DELTA-ERD = erd_slow - erd_fast
%   p5_flag    : 1 if pre-stim power in either sub-band < 5th percentile (session)
%
% featGA has the same field names but averaged across channels (omitnan),
% returned as [nTrials x 1] column vectors. p5_flag in featGA is a logical
% OR across channels (any flagged channel flags the trial).
%
% Inputs:
%   f        : [1 x nFreq]  frequency vector (Hz)
%   prePxx   : [nChan x nFreq x nTrials]  pre-stimulus Welch PSD
%   postPxx  : [nChan x nFreq x nTrials]  post-stimulus Welch PSD
%   alpha    : struct with fields  alpha_hz [1x2], slow_hz [1x2], fast_hz [1x2]
%              alpha_hz should be [8 12]; CoG is computed over this same band.
%   logf     : (optional) MATLAB file handle for spec_logmsg
%
% Notes:
%   - CoG is computed over alpha_hz = [8, 12] Hz (same as slow+fast union).
%     The fixed 10 Hz boundary in psi_cog reflects the slow/fast split point,
%     not an assumption about the CoG range.
%   - prePxx and postPxx are expected to arrive pre-windowed by the caller:
%       pre  : [-1000, -100] ms relative to stimulus onset
%       post : [+100,  +800] ms relative to stimulus onset
%   - TVI_alpha (between-subjects) and phase (Hilbert-based) are computed
%     upstream and are not part of this function.

if nargin < 4
    error('spec_compute_interaction_metrics:MissingArg', ...
        'spec_compute_interaction_metrics requires at least 4 arguments (f, prePxx, postPxx, alpha); got %d.', nargin);
end

if nargin < 5, logf = 1; end

if isempty(prePxx) || isempty(postPxx)
    error('spec_compute_interaction_metrics:EmptyPxx', ...
        'prePxx/postPxx must not be empty (isempty(prePxx) = %d, isempty(postPxx) = %d).', isempty(prePxx), isempty(postPxx));
end

if ~isequal(size(prePxx), size(postPxx))
    error('spec_compute_interaction_metrics:SizeMismatch', ...
        'prePxx and postPxx must have identical shape. Got prePxx = %s, postPxx = %s.', mat2str(size(prePxx)), mat2str(size(postPxx)));
end

if numel(f) ~= size(prePxx, 2)
    error('spec_compute_interaction_metrics:FreqMismatch', ...
        'numel(f) = %d does not match size(prePxx, 2) = %d.', numel(f), size(prePxx, 2));
end

slow_hz  = alpha.slow_hz;   % [8  10]
fast_hz  = alpha.fast_hz;   % [10 12]
alpha_hz = alpha.alpha_hz;  % [8  12]

f = f(:)';   % ensure row vector for indexing

idxS = (f >= slow_hz(1))  & (f <= slow_hz(2));
idxF = (f >= fast_hz(1))  & (f <= fast_hz(2));
idxA = (f >= alpha_hz(1)) & (f <= alpha_hz(2));

if ~any(idxS) || ~any(idxF) || ~any(idxA)
    error('spec_compute_interaction_metrics:BadBands', ...
        'No frequency bins found in one or more alpha sub-bands. Check alpha_hz / slow_hz / fast_hz vs f range [%.1f %.1f].', ...
        f(1), f(end));
end

fS = f(idxS);
fF = f(idxF);
fA = f(idxA);

nChan = size(prePxx, 1);
nTr   = size(prePxx, 3);

% ---------------------------------------------------------------
% Band power via trapz  ->  [nChan x nTrials]
% ensure_2d guards against squeeze dropping a dimension when
% nChan=1 or nTrials=1, using explicit target shape [nChan, nTr].
% ---------------------------------------------------------------
pow_pre_slow = band_power(fS, prePxx, idxS, nChan, nTr, 'pow_pre_slow');
pow_pre_fast = band_power(fF, prePxx, idxF, nChan, nTr, 'pow_pre_fast');
pow_pre_alpha = band_power(fA, prePxx, idxA, nChan, nTr, 'pow_pre_alpha');
pow_post_slow = band_power(fS, postPxx, idxS, nChan, nTr, 'pow_post_slow');
pow_post_fast = band_power(fF, postPxx, idxF, nChan, nTr, 'pow_post_fast');
pow_post_alpha = band_power(fA, postPxx, idxA, nChan, nTr, 'pow_post_alpha');

% ---------------------------------------------------------------
% Noise floor epsilon_0  (see compute_noise_floor below)
% ---------------------------------------------------------------
eps0 = compute_noise_floor(f, prePxx);
spec_logmsg(logf, '[INTERACT] eps0 (noise floor) = %.4g uV^2/Hz', eps0);

% ---------------------------------------------------------------
% ERD family (signed; negative = power decrease = desynchronization)
% Each sub-band / total-alpha normalizes by its own pre-stim baseline
% ---------------------------------------------------------------
erd_slow = (pow_post_slow - pow_pre_slow) ./ (pow_pre_slow + eps0);
erd_fast = (pow_post_fast - pow_pre_fast) ./ (pow_pre_fast + eps0);
erd_pow_alpha_total = (pow_post_alpha - pow_pre_alpha) ./ (pow_pre_alpha + eps0);
delta_erd = erd_slow - erd_fast;

% ---------------------------------------------------------------
% p5_flag: mark trials where pre-stim power in either sub-band
% falls below the 5th percentile of the session distribution.
% Computed across all valid (non-NaN) [chan x trial] values.
% ---------------------------------------------------------------
all_slow = pow_pre_slow(~isnan(pow_pre_slow));
all_fast = pow_pre_fast(~isnan(pow_pre_fast));

if isempty(all_slow) || isempty(all_fast)
    p5_flag = zeros(nChan, nTr);
    spec_logmsg(logf, '[INTERACT][WARN] All pre-stim ower values are NaN; p5_flag set to zero.');
else
    thr_slow = prctile(all_slow, 5);
    thr_fast = prctile(all_fast, 5);
    p5_flag = double((pow_pre_slow < thr_slow) | (pow_pre_fast < thr_fast));
    spec_logmsg(logf, '[INTERACT] p5 thresholds: slow = %.4g fast = %.4g flagged %d / %d cells', ...
        thr_slow, thr_fast, sum(p5_flag(:)), numel(p5_flag));
end

% ---------------------------------------------------------------
% Pack per-channel x trial struct
% ---------------------------------------------------------------
featChan = struct( ...
    'erd_slow', erd_slow, ...
    'erd_fast', erd_fast, ...
    'erd_pow_alpha_total', erd_pow_alpha_total, ...
    'delta_erd', delta_erd, ...
    'p5_flag', p5_flag ...
);

% ---------------------------------------------------------------
% GA: mean across channels (omitnan), normalized to [nTrials x 1].
% Exception: p5_flag uses logical OR across channels -- a trial is
% flagged if ANY channel meets the threshold, not on the average.
% ---------------------------------------------------------------
fn = fieldnames(featChan);
featGA = struct();
for i = 1:numel(fn)
    if strcmp(fn{i}, 'p5_flag')
        featGA.p5_flag = double(any(featChan.p5_flag, 1))'; % [nTrials x 1]
    else
        m = mean(featChan.(fn{i}), 1, 'omitnan'); % [1 x nTrials]
        featGA.(fn{i}) = m(:); % [nTrials x 1]
    end
end
end

% ================================================================
% Local: noise-floor epsilon_0
% Try the 45-55 Hz "quiet band" first; fall back to a small
% fraction of the median total power if the band in unavailable
% ================================================================
function eps0 = compute_noise_floor(f, Pxx)
lo = 45; hi = 55;
idxQ = (f >= lo) & (f <= hi);

if sum(idxQ) >= 2
    qvals = Pxx(:, idxQ, :);
    eps0 = median(qvals(~isnan(qvals)), 'all');
    if isnan(eps0) || eps0 <= 0
        eps0 = median(Pxx(~isnan(Pxx)), 'all') * 1e-3;
    end
else
    eps0 = median(Pxx(~isnan(Pxx)), 'all') * 1e-3;
end

eps0 = max(eps0, 1e-12); % hard lower bound
end

% ================================================================
% Local: band power via trapz, shape-checked with a caller-supplied
% name so a failure identifies exactly which (band x window) broke.
% ================================================================
function X = band_power(fBand, Pxx, idxBand, nChan, nTr, label)
try
    X = squeeze(trapz(fBand, Pxx(:, idxBand, :), 2));
    if ~isequal(size(X), [nChan, nTr])
        X = reshape(X, nChan, nTr);
    end
catch ME
    error('spec_compute_interaction_metrics:BandPowerFailed', ...
        '%s failed: Pxx size = %s, idxBand true-count = %d, target shape = [%d %d]. Underlying error: %s', ...
        label, mat2str(size(Pxx)), sum(idxBand), nChan, nTr, ME.message);
end
end