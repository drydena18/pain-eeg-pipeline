function tvi = spec_compute_tvi_alpha(bi_pre, logf)
% SPEC_COMPUTE_TVI_ALPHA  Temporal variability index of a pre_sf_balance sequence
% V 1.1.0
%
% Python twin: src_prestim.src_compute_tvi_alpha.
%
%   TVI_alpha = MSSD / Var
%   MSSD = mean((b(k+1) - b(k)).^2)
%   Var  = unbiased sample variance, var(b, 0)
%
% Range [0, 4]: ~0 slowly drifting state, ~2 trial-to-trial independence,
% > 2 alternation between successive trials.
%
% V1.1.0: variance floor is 1e-12 (was eps) to match Python; pass logf = -1
% to suppress logging (used for the per-channel loop); docstring range
% corrected from [0, 2].
%
% Output struct: tvi_alpha, mssd_bi_pre, var_bi_pre, n_trials_bi

if nargin < 2, logf = 1; end
quiet = isnumeric(logf) && isscalar(logf) && logf < 0;

b    = bi_pre(:);
b_ok = b(~isnan(b));
K    = numel(b_ok);

tvi = struct('tvi_alpha', nan, 'mssd_bi_pre', nan, 'var_bi_pre', nan, 'n_trials_bi', K);
if K < 3
    if ~quiet
        spec_logmsg(logf, '[TVI] Only %d valid trials (need >= 3); TVI_alpha = NaN.', K);
    end
    return;
end

mssd = mean(diff(b_ok) .^ 2);
vr   = var(b_ok, 0);
tvi.mssd_bi_pre = mssd;
tvi.var_bi_pre  = vr;
if vr >= 1e-12
    tvi.tvi_alpha = mssd / vr;
end
if ~quiet
    spec_logmsg(logf, '[TVI] n=%d  TVI_alpha=%.4f  MSSD=%.4g  Var=%.4g', K, tvi.tvi_alpha, mssd, vr);
end
end
