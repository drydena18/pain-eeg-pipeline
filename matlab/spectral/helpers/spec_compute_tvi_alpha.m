function tvi = spec_compute_tvi_alpha(bi_pre, logf)
% SPEC_COMPUTE_TVI_ALPHA  Temporal variability index from the GA BI_pre sequence
% V 1.0.0
%
% TVI_alpha quantifies how dynamically a subject's pre-stimulus alpha state
% fluctuates across the session, independently of total amplitude variance.
%
%   TVI_alpha = MSSD / Var(BI_pre)
%
% Where:
%   MSSD = (1 / K-1) * sum( (BI_pre(k+1) - BI_pre(k))^2 )
%            mean square successive difference — sensitive to the
%            trial-by-trial autocorrelation structure of the sequence
%   Var  = unbiased sample variance of the BI_pre sequence
%            normalises MSSD so the result is independent of total swing
%
% Range [0, 2]:
%   Near 0  : slowly varying / rigid pre-stimulus state across the session
%   Near 2  : rapidly alternating from trial to trial (theoretical maximum
%              for a perfectly alternating +1/-1 sequence)
%
% Requires >= 3 valid (non-NaN) trials; returns NaN for all scalar fields
% otherwise.
%
% Inputs:
%   bi_pre : [nTrials x 1] or [1 x nTrials]  GA BI_pre sequence
%   logf   : (optional) MATLAB log file handle
%
% Output:
%   tvi : struct
%           .tvi_alpha    scalar  NaN if < 3 valid trials or Var(bi_pre) ~ 0
%           .mssd_bi_pre  scalar  NaN if < 3 valid trials
%           .var_bi_pre   scalar  NaN if < 3 valid trials
%           .n_trials_bi  integer count of valid (non-NaN) trials used

if nargin < 2, logf = 1; end

b    = bi_pre(:);
b_ok = b(~isnan(b));
K    = numel(b_ok);

tvi = struct( ...
    'tvi_alpha',   nan, ...
    'mssd_bi_pre', nan, ...
    'var_bi_pre',  nan, ...
    'n_trials_bi', K   ...
);

if K < 3
    spec_logmsg(logf, '[TVI] Only %d valid BI_pre trials (need >= 3); TVI_alpha = NaN.', K);
    return;
end

diffs = diff(b_ok);
mssd  = mean(diffs .^ 2);
vr    = var(b_ok, 0);   % unbiased; denominator K-1

tvi.mssd_bi_pre = mssd;
tvi.var_bi_pre  = vr;

if vr > eps
    tvi.tvi_alpha = mssd / vr;
else
    spec_logmsg(logf, '[TVI] Var(BI_pre) ~ 0 for this subject; TVI_alpha = NaN.');
end

spec_logmsg(logf, '[TVI] n=%d  TVI_alpha=%.4f  MSSD=%.4g  Var=%.4g', ...
    K, tvi.tvi_alpha, tvi.mssd_bi_pre, tvi.var_bi_pre);
end