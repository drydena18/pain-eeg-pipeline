function q = spec_percentile_linear(v, p)
% SPEC_PERCENTILE_LINEAR  NaN-aware percentile, numpy 'linear' definition
% V 1.0.0
%
% Python twin: numpy.nanpercentile(v, p) (default method 'linear'):
%   h = (n - 1) * p / 100 on the sorted values, linear interpolation.
% MATLAB's prctile uses a different interpolation rule, so it is not used.

v = sort(double(v(~isnan(v(:)))));
n = numel(v);
if n == 0
    q = NaN;
    return;
end
h  = (n - 1) * p / 100;
lo = floor(h);
hi = min(lo + 1, n - 1);
q  = v(lo + 1) + (h - lo) * (v(hi + 1) - v(lo + 1));
end
