function Y = spec_zero_phase_filter(X, h)
% SPEC_ZERO_PHASE_FILTER  Zero-phase FIR along dim 2 of [nUnit x nTime x nTr]
% V 1.0.0
%
% Python twin: src_spectral.src_zero_phase_filter.
%   - mirror-pad (edge sample not repeated) by min(N, nTime - 1) samples
%   - FFT convolution with the symmetric odd-length kernel, centred ('same')
%   - crop back to nTime
%
% Processed in column chunks to bound memory.

h = double(h(:));
N = numel(h);
if mod(N, 2) == 0
    error('spec_zero_phase_filter:EvenKernel', 'Kernel length must be odd.');
end

sz = size(X);
if numel(sz) < 3, sz(3) = 1; end
nT = sz(2);
pad = min(N, nT - 1);
M = (N - 1) / 2;

C = reshape(permute(double(X), [2 1 3]), nT, []);   % [nT x nUnit*nTr]
nCol = size(C, 2);
L = nT + 2 * pad;
nfft = 2 ^ nextpow2(L + N - 1);
H = fft(h, nfft);

Yc = zeros(nT, nCol);
chunk = 512;
for c0 = 1:chunk:nCol
    cols = c0:min(c0 + chunk - 1, nCol);
    Cp = [C(pad+1:-1:2, cols); C(:, cols); C(end-1:-1:end-pad, cols)];
    Yf = real(ifft(fft(Cp, nfft) .* H));
    Ys = Yf(M+1:M+L, :);                 % 'same' (centred) part
    Yc(:, cols) = Ys(pad+1:pad+nT, :);
end

Y = permute(reshape(Yc, nT, sz(1), sz(3)), [2 1 3]);
end
