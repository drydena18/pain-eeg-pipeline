function r = pp_eig_rank(data, tol)
%PP_EIG_RANK Rank estimate from covariance eigenvalues, computed in chunks
%   so long continuous recordings don't need a full double copy in memory.
%   Used only as a cross-check on the bookkeeping rank.
if nargin < 2, tol = 1e-7; end
n = size(data, 1);
data = reshape(data, n, []);
T = size(data, 2);
mu = mean(data, 2, 'double');
C = zeros(n);
step = 1e5;
for s = 1:step:T
    e = min(s + step - 1, T);
    X = double(data(:, s:e)) - mu;
    C = C + X * X';
end
ev = eig((C + C') / 2);
r = sum(ev > tol * max(ev));
end