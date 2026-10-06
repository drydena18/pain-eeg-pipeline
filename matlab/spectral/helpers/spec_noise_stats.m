function st = spec_noise_stats(preSlow, preFast, eps0Frac, pPct)
% SPEC_NOISE_STATS  Per-unit ERD guard eps0 and p5_flag thresholds
% V 1.0.0
%
% Python twin: src_erd.src_compute_noise_stats. Unit = channel here, ROI in
% source space; statistics are pooled across trials WITHIN each unit:
%
%   eps0     = max(eps0Frac * median([preSlow, preFast]), 1e-12)
%   thr_slow = pPct-th percentile of preSlow   (numpy 'linear')
%   thr_fast = pPct-th percentile of preFast
%
% Replaces compute_noise_floor in spec_compute_interaction_metrics.m, whose
% 45-55 Hz quiet band lies in the stopband of the 40 Hz preprocessing
% low-pass and whose statistics were pooled across all channels.
%
% Inputs:
%   preSlow, preFast : [nUnit x nTr] pre-stim band power
%   eps0Frac         : default 1e-3
%   pPct             : default 5
%
% Output:
%   st.eps0, st.thr_slow, st.thr_fast : [nUnit x 1]

if nargin < 3 || isempty(eps0Frac), eps0Frac = 1e-3; end
if nargin < 4 || isempty(pPct),     pPct = 5;        end

nU = size(preSlow, 1);
st = struct('eps0', nan(nU, 1), 'thr_slow', nan(nU, 1), 'thr_fast', nan(nU, 1));
for u = 1:nU
    s = preSlow(u, :);
    q = preFast(u, :);
    st.thr_slow(u) = spec_percentile_linear(s, pPct);
    st.thr_fast(u) = spec_percentile_linear(q, pPct);
    pooled = [s(:); q(:)];
    pooled = pooled(~isnan(pooled));
    if isempty(pooled)
        e0 = 1e-12;
    else
        e0 = median(pooled) * eps0Frac;
    end
    st.eps0(u) = max(e0, 1e-12);
end
end
