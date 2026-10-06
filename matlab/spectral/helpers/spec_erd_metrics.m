function erd = spec_erd_metrics(pre, post, st, useP5)
% SPEC_ERD_METRICS  ERD family (fractional pre -> post change) + p5_flag
% V 1.1.0
%
% Python twin: src_erd.src_compute_erd_metrics.
%   erd_b     = (post_b - pre_b) ./ (pre_b + eps0),  b in {slow, fast, alpha}
%   erd_asym  = erd_slow - erd_fast   (was delta_erd; delta_ = post - pre only)
%   p5_flag   = pre_slow < thr_slow | pre_fast < thr_fast   (if useP5)
%
% Inputs:
%   pre, post : structs with .slow .fast .alpha  [nUnit x nTr] (or [nUnit x 1])
%   st        : from spec_noise_stats (per-unit [nUnit x 1] vectors)
%   useP5     : false for subject-level (trial-mean) rows
%
% Output:
%   erd : struct erd_slow, erd_fast, erd_pow_alpha_total, erd_asym [, p5_flag]

if nargin < 4, useP5 = true; end

e0 = st.eps0(:);
erd = struct();
erd.erd_slow            = (post.slow  - pre.slow)  ./ (pre.slow  + e0);
erd.erd_fast            = (post.fast  - pre.fast)  ./ (pre.fast  + e0);
erd.erd_pow_alpha_total = (post.alpha - pre.alpha) ./ (pre.alpha + e0);
erd.erd_asym            = erd.erd_slow - erd.erd_fast;

if useP5
    erd.p5_flag = double((pre.slow < st.thr_slow(:)) | (pre.fast < st.thr_fast(:)));
end
end
