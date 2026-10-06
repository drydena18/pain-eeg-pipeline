function cog = spec_cog_from_psd(f, P, band)
% SPEC_COG_FROM_PSD  Spectral centre of gravity over band, along dim 2
% V 1.0.0
%
% Python twin: src_spectral.src_cog_from_psd.
%   cog = trapz(f, P .* f) / trapz(f, P) over bins inside band
% NaN when fewer than 2 bins fall in the band or band power <= 0.
%
% Inputs:
%   f    : [1 x nF]
%   P    : [nUnit x nF x nTr] (or [nUnit x nF])
%   band : [lo hi] Hz
%
% Output:
%   cog  : [nUnit x nTr]

f = f(:)';
sz = size(P);
if numel(sz) < 3, sz(3) = 1; end
idx = f >= band(1) - 1e-9 & f <= band(2) + 1e-9;
if nnz(idx) < 2
    cog = nan(sz(1), sz(3));
    return;
end
fb = f(idx);
Pb = P(:, idx, :);
den = trapz(fb, Pb, 2);
num = trapz(fb, Pb .* reshape(fb, 1, [], 1), 2);
cog = reshape(num ./ den, sz(1), sz(3));
den = reshape(den, sz(1), sz(3));
cog(~(isfinite(den) & den > 0)) = NaN;
end
