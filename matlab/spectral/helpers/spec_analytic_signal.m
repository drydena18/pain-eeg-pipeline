function Z = spec_analytic_signal(X)
% SPEC_ANALYTIC_SIGNAL  FFT-based analytic signal along dim 2
% V 1.0.0
%
% Python twin: src_spectral.src_analytic_signal. Same definition as
% MATLAB's hilbert() / scipy.signal.hilbert, written out so it needs no
% toolbox and works along dim 2 of [nUnit x nTime x nTr].

n = size(X, 2);
g = zeros(1, n);
g(1) = 1;
if mod(n, 2) == 0
    g(2:n/2) = 2;
    g(n/2 + 1) = 1;
else
    g(2:(n + 1)/2) = 2;
end
Z = ifft(fft(double(X), [], 2) .* g, [], 2);
end
