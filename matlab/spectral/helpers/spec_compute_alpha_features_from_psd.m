function feat = spec_compute_alpha_features_from_psd(f, Pxx, alpha)
% f:    [1 x nFreq]
% Pxx:  [nChan x nFreq x nTrials]
%
% Returns fields [nChan x nTrials] (or [1 x nTrials] for GA wrapper)

eps0 = 1e-12;

slow = alpha.slow_hz;
fast = alpha.fast_hz;
abnd = alpha.alpha_hz;

idxS = (f >= slow(1)) & (f <= slow(2));
idxF = (f >= fast(1)) & (f <= fast(2));
idxA = (f >= abnd(1)) & (f <= abnd(2));

% bandpower via trapz
pow_slow = squeeze(trapz(f(idxS), Pxx(:, idxS, :), 2)); % [chan x trial]
pow_fast = squeeze(trapz(f(idxF), Pxx(:, idxF, :), 2));
pow_alpha = squeeze(trapz(f(idxA), Pxx(:, idxA, :), 2));

% PAF CoG within alpha band
nums = squeeze(trapz(f(idxA), Pxx(:, idxA, :) .* reshape(f(idxA), [1 sum(idxA) 1]), 2));
den = pow_alpha + eps0;
paf_cog = num ./ den;

% Interaction family
sf_ratio    = pow_slow ./ (pow_fast + eps0);
sf_logratio = log(pow_slow + eps0) - log(pow_fast + eps0);
sf_balance  = (pow_slow - pow_fast) ./ (pow_slow + pow_fast + eps0);
slow_frac   = pow_slow ./ (pow_slow + pow_fast + eps0);

% Relative within alpha
rel_slow = pow_slow ./ (pow_alpha + eps0);
rel_fast = pow_fast ./ (pow_alpha + eps0);

feat = struct();
feat.paf_cog_hz = paf_cog;
feat.pow_slow_alpha = pow_slow;
feat.pow_fast_alpha = pow_fast;
feat.pow_alpha_total = pow_alpha;
feat.rel_slow_alpha = rel_slow;
feat.rel_fast_alpha = rel_fast;
feat.sf_ratio = sf_ratio;
feat.sf_logratio = sf_logratio;
feat.sf_balance = sf_balance;
feat.slow_alpha_frac = slow_frac;
end