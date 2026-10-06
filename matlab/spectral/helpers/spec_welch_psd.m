function [f, P] = spec_welch_psd(X, fs, fmin, fmax, windowSec, overlapFrac, dfTarget)
% SPEC_WELCH_PSD  Shared Welch PSD along dim 2 of [nUnit x nSamp x nTr]
% V 1.0.0
%
% Python twin: src_spectral.src_psd_welch (scipy.signal.welch with explicit
% settings). Written out so segmentation, window and scaling are identical:
%
%   segment  = min(floor(windowSec * fs), nSamp) samples (>= 8)
%   step     = max(1, floor(segment * (1 - overlapFrac)))
%   nfft     = max(2^ceil(log2(segment)), ceil(fs / dfTarget))
%   window   = PERIODIC Hamming (scipy 'hamming'); pwelch would use symmetric
%   no detrending, density scaling 1 / (fs * sum(w.^2)), one-sided,
%   mean across segments
%
% Outputs:
%   f : [1 x nF] Hz within [fmin, fmax]
%   P : [nUnit x nF x nTr]

if nargin < 5 || isempty(windowSec),   windowSec = 2.0;   end
if nargin < 6 || isempty(overlapFrac), overlapFrac = 0.5; end
if nargin < 7 || isempty(dfTarget),    dfTarget = 0.25;   end

sz = size(X);
if numel(sz) < 3, sz(3) = 1; end
nS = sz(2);

nwin  = max(min(floor(windowSec * fs + 1e-9), nS), 8);
nstep = max(1, floor(nwin * (1 - overlapFrac) + 1e-9));
nfft  = max(2 ^ ceil(log2(nwin)), ceil(fs / dfTarget - 1e-9));

w = 0.54 - 0.46 * cos(2 * pi * (0:nwin-1)' / nwin);   % periodic Hamming
scale = 1 / (fs * sum(w .^ 2));
starts = 1:nstep:(nS - nwin + 1);
nOne = floor(nfft / 2) + 1;

C = reshape(permute(double(X), [2 1 3]), nS, []);   % [nS x nUnit*nTr]
nCol = size(C, 2);
Pc = zeros(nOne, nCol);

chunk = 512;
for c0 = 1:chunk:nCol
    cols = c0:min(c0 + chunk - 1, nCol);
    acc = zeros(nOne, numel(cols));
    for s = starts
        F = fft(C(s:s+nwin-1, cols) .* w, nfft);
        acc = acc + abs(F(1:nOne, :)) .^ 2;
    end
    Pc(:, cols) = acc * scale / numel(starts);
end

if mod(nfft, 2) == 0
    Pc(2:end-1, :) = 2 * Pc(2:end-1, :);
else
    Pc(2:end, :) = 2 * Pc(2:end, :);
end

fAll = (0:nOne-1)' * fs / nfft;
keep = fAll >= fmin - 1e-9 & fAll <= fmax + 1e-9;
f = fAll(keep)';
P = permute(reshape(Pc(keep, :), nnz(keep), sz(1), sz(3)), [2 1 3]);
end
