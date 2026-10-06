function PT = spec_band_power_tcs(X, fs, bands, transBw)
% SPEC_BAND_POWER_TCS  Filter-Hilbert band-power time courses, full epoch
% V 1.0.0
%
% Python twin: src_alpha_features.src_compute_band_power_tcs /
% src_spectral.src_band_power_tc.
%
%   power(t) = |analytic(x_band)(t)|^2 / 2
%
% The window mean of power(t) is the variance of the band-limited signal,
% i.e. the one-sided PSD integrated over the band (Parseval). Filter the
% FULL epoch and crop afterwards so edge transients stay outside windows.
%
% Inputs:
%   X       : [nUnit x nTime x nTr]
%   fs      : Hz
%   bands   : struct, each field a [lo hi] band in Hz (e.g. slow/fast/alpha)
%   transBw : FIR transition bandwidth in Hz
%
% Output:
%   PT : struct with the same field names, each [nUnit x nTime x nTr]

names = fieldnames(bands);
PT = struct();
for i = 1:numel(names)
    b = bands.(names{i});
    h = spec_fir_bandpass(fs, b(1), b(2), transBw);
    if numel(h) > size(X, 2)
        error('spec_band_power_tcs:FilterTooLong', ...
            'FIR length %d exceeds epoch length %d; raise features.band_power.trans_bw_hz.', ...
            numel(h), size(X, 2));
    end
    Z = spec_analytic_signal(spec_zero_phase_filter(X, h));
    PT.(names{i}) = abs(Z) .^ 2 / 2;
end
end
