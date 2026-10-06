function h = spec_fir_bandpass(fs, lo, hi, transBw)
% SPEC_FIR_BANDPASS  Shared linear-phase windowed-sinc band-pass FIR
% V 1.0.0
%
% Python twin: src_spectral.src_fir_bandpass. Written out explicitly (no
% fir1 / pop_eegfiltnew) so MATLAB and Python use the identical kernel:
%
%   N = floor(3.3 * fs / transBw + 0.5), forced odd   (Hamming rule)
%   h = (lowpass(hi) - lowpass(lo)) .* symmetric Hamming(N)
%   normalised to unit gain at (lo + hi) / 2
%
% The -6 dB points sit exactly on lo and hi, so 8-10 and 10-12 Hz meet at
% 10 Hz. transBw is centred on each edge.
%
% Output:
%   h : [1 x N] taps

if ~(lo > 0 && lo < hi && hi < fs / 2)
    error('spec_fir_bandpass:BadBand', 'Bad band [%g %g] Hz for fs = %g Hz.', lo, hi, fs);
end

N = floor(3.3 * fs / transBw + 0.5);
if mod(N, 2) == 0
    N = N + 1;
end
n = 0:N-1;
k = n - (N - 1) / 2;

win = 0.54 - 0.46 * cos(2 * pi * n / (N - 1));
h = (local_lowpass(hi, fs, k) - local_lowpass(lo, fs, k)) .* win;

fc = (lo + hi) / 2;
g = abs(sum(h .* exp(-2i * pi * fc * n / fs)));
h = h / g;
end

function y = local_lowpass(fc, fs, k)
w = 2 * fc / fs;
x = w * k;
s = ones(size(x));
nz = x ~= 0;
s(nz) = sin(pi * x(nz)) ./ (pi * x(nz));
y = w * s;
end
