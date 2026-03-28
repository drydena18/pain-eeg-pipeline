function [featChan, featGA] = spec_compute_interaction_metrics(f, prePxx, postPxx, alpha, logf)
% SPEC_COMPUTE_INTERACTION_METRICS  Pre-stimulus alpha interaction metrics
% V 1.0.0
%
% Computes the following per-channel x trial fields:
%   bi_pre     : Balance Index [-1,+1]  (slow-fast)/(slow+fast+eps0)
%   lr_pre     : Log Ratio (unbounded)  log(slow+eps0) - log(fast+eps0)
%   cog_pre    : Centre of Gravity [Hz] over full alpha band [alpha_hz]
%   psi_cog    : Interaction term  bi_pre x (cog_pre - 10)
%   erd_slow   : Slow-alpha ERD  (post-pre)/pre  (negative = desynchronisation)
%   erd_fast   : Fast-alpha ERD  (post-pre)/pre
%   delta_erd  : DELTA-ERD = erd_slow - erd_fast
%   p5_flag    : 1 if pre-stim power in either sub-band < 5th percentile (session)
%
% featGA has the same field names but averaged across channels (omitnan),
% returned as [nTrials x 1] column vectors.
%
% Inputs:
%   f        : [1 x nFreq]  frequency vector (Hz)
%   prePxx   : [nChan x nFreq x nTrials]  pre-stimulus Welch PSD
%   postPxx  : [nChan x nFreq x nTrials]  post-stimulus Welch PSD
%   alpha    : struct with fields  alpha_hz [1x2], slow_hz [1x2], fast_hz [1x2]
%   logf     : (optional) MATLAB file handle for spec_logmsg

if nargin < 5, logf = 1; end

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

% ---------------------------------------------------------------
% Band power via trapz  ->  [nChan x nTrials]
% Note: squeeze collapses the frequency dimension (dim 2) after trapz.
% When nChan=1 or nTrials=1, squeeze may reduce rank; guard below.
% ---------------------------------------------------------------
pow_pre_slow  = ensure_2d(squeeze(trapz(fS, prePxx(:, idxS, :), 2)));
pow_pre_fast  = ensure_2d(squeeze(trapz(fF, prePxx(:, idxF, :), 2)));
pow_post_slow = ensure_2d(squeeze(trapz(fS, postPxx(:, idxS, :), 2)));
pow_post_fast = ensure_2d(squeeze(trapz(fF, postPxx(:, idxF, :), 2)));
pow_pre_alpha = ensure_2d(squeeze(trapz(fA, prePxx(:, idxA, :), 2)));

nChan = size(prePxx, 1);
nTr   = size(prePxx, 3);

% Verify shapes after squeeze
expected = [nChan, nTr];
for chk = {pow_pre_slow, pow_pre_fast, pow_post_slow, pow_post_fast, pow_pre_alpha}
    if ~isequal(size(chk{1}), expected)
        error('spec_compute_interaction_metrics:ShapeMismatch', ...
            'Band-power array shape %s != expected [%d %d]. Check prePxx/postPxx dimensions.', ...
            mat2str(size(chk{1})), expected(1), expected(2));
    end
end

% ---------------------------------------------------------------
% Noise floor epsilon_0  (see compute_noise_floor below)
% ---------------------------------------------------------------
eps0 = compute_noise_floor(f, prePxx);
spec_logmsg(logf, '[INTERACT] eps0 (noise floor) = %.4g uV^2/Hz', eps0);

% ---------------------------------------------------------------
% BI_pre and LR_pre
% ---------------------------------------------------------------
bi_pre = (pow_pre_slow - pow_pre_fast) ./ (pow_pre_slow + pow_pre_fast + eps0);
lr_pre = log(pow_pre_slow + eps0) - log(pow_pre_fast + eps0);

% ---------------------------------------------------------------
% CoG_pre  (frequency-weighted centroid over full alpha band)
% Uses trapz to match the band-power integration above
% ---------------------------------------------------------------
fA_3d   = reshape(fA, [1, numel(fA), 1]);   % [1 x nAlphaBins x 1] for broadcasting
num_cog = ensure_2d(squeeze(trapz(fA, bsxfun(@times, prePxx(:, idxA, :), fA_3d), 2)));
cog_pre = num_cog ./ (pow_pre_alpha + eps0);

% Interaction term: BI_pre x (CoG_pre - 10 Hz boundary)
psi_cog = bi_pre .* (cog_pre - 10);

% ---------------------------------------------------------------
% ERD (signed; negative = power decrease = desynchronisation)
%   ERD_slow = (post_slow - pre_slow) / pre_slow
%   ERD_fast = (post_fast - pre_fast) / pre_fast
% ---------------------------------------------------------------
erd_slow  = (pow_post_slow - pow_pre_slow) ./ (pow_pre_slow + eps0);
erd_fast  = (pow_post_fast - pow_pre_fast) ./ (pow_pre_fast + eps0);
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
    spec_logmsg(logf, '[INTERACT][WARN] All pre-stim power values are NaN; p5_flag set to zero.');
else
    thr_slow = prctile(all_slow, 5);
    thr_fast = prctile(all_fast, 5);
    p5_flag  = double((pow_pre_slow < thr_slow) | (pow_pre_fast < thr_fast));
    spec_logmsg(logf, '[INTERACT] p5 thresholds: slow=%.4g fast=%.4g  flagged %d / %d cells', ...
        thr_slow, thr_fast, sum(p5_flag(:)), numel(p5_flag));
end

% ---------------------------------------------------------------
% Pack per-channel x trial struct
% ---------------------------------------------------------------
featChan = struct( ...
    'bi_pre',    bi_pre,    ...
    'lr_pre',    lr_pre,    ...
    'cog_pre',   cog_pre,   ...
    'psi_cog',   psi_cog,   ...
    'erd_slow',  erd_slow,  ...
    'erd_fast',  erd_fast,  ...
    'delta_erd', delta_erd, ...
    'p5_flag',   p5_flag    ...
);

% ---------------------------------------------------------------
% GA: mean across channels (omitnan), normalised to [nTrials x 1]
% ---------------------------------------------------------------
fn = fieldnames(featChan);
featGA = struct();
for i = 1:numel(fn)
    m = mean(featChan.(fn{i}), 1, 'omitnan');  % [1 x nTrials]
    featGA.(fn{i}) = m(:);                      % [nTrials x 1]
end
end

% ================================================================
%  Local: noise-floor epsilon_0
%  Try the 45-55 Hz "quiet band" first; fall back to a small
%  fraction of the median total power if the band is unavailable.
% ================================================================
function eps0 = compute_noise_floor(f, Pxx)
lo = 45; hi = 55;
idxQ = (f >= lo) & (f <= hi);

if sum(idxQ) >= 2
    qvals = Pxx(:, idxQ, :);
    eps0  = median(qvals(~isnan(qvals)), 'all');
    if isnan(eps0) || eps0 <= 0
        eps0 = median(Pxx(~isnan(Pxx)), 'all') * 1e-3;
    end
else
    eps0 = median(Pxx(~isnan(Pxx)), 'all') * 1e-3;
end

eps0 = max(eps0, 1e-12);   % hard lower bound
end

% ================================================================
%  Local: ensure a band-power result is [nChan x nTrials]
%  squeeze() on a 1-channel or 1-trial array may drop a dimension.
% ================================================================
function X = ensure_2d(X)
if isvector(X)
    % cannot infer correct shape without context; leave as-is and let
    % the shape-check in the caller catch genuine mismatches
end
if isscalar(X)
    X = reshape(X, 1, 1);
end
end