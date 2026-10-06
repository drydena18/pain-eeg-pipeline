function [ph, itc] = spec_slow_phase_features(X, timesSec, fs, slowBand, transBw, postRef, postWin)
% SPEC_SLOW_PHASE_FEATURES  Slow-alpha Hilbert phase at t = 0 / postRef + ITC
% V 1.0.0
%
% Python twins: src_prestim.src_compute_prestim_phase,
% src_poststim.src_compute_poststim_phase, src_poststim.src_compute_itc.
% Same FIR as band power (spec_fir_bandpass), applied to the FULL epoch.
%
% Inputs:
%   X        : [nUnit x nTime x nTr] (pass mean(X, 3) for trial-mean phase)
%   timesSec : [1 x nTime]
%   fs, slowBand [lo hi], transBw, postRef (s), postWin [tmin tmax] (s)
%
% Outputs:
%   ph  : struct slow_phase, sin_phase, cos_phase,
%         slow_phase_post, sin_phase_post, cos_phase_post   [nUnit x nTr]
%   itc : struct itc_mean, itc_peak, itc_peak_latency_ms     [nUnit x 1]
%         (inter-trial phase coherence over postWin; first maximum)

timesSec = double(timesSec(:))';
nU = size(X, 1);
nTr = size(X, 3);

h = spec_fir_bandpass(fs, slowBand(1), slowBand(2), transBw);
Z = spec_analytic_signal(spec_zero_phase_filter(X, h));

[~, i0] = min(abs(timesSec - 0));
[~, iR] = min(abs(timesSec - postRef));

ph = struct();
p0 = reshape(angle(Z(:, i0, :)), nU, nTr);
pR = reshape(angle(Z(:, iR, :)), nU, nTr);
ph.slow_phase      = p0;
ph.sin_phase       = sin(p0);
ph.cos_phase       = cos(p0);
ph.slow_phase_post = pR;
ph.sin_phase_post  = sin(pR);
ph.cos_phase_post  = cos(pR);

if nargout > 1
    m = spec_window_mask(timesSec, postWin(1), postWin(2), fs);
    tPost = timesSec(m);
    U = exp(1i * angle(Z(:, m, :)));
    I = abs(mean(U, 3));                           % [nUnit x nPost]
    [pk, ip] = max(I, [], 2);
    itc = struct();
    itc.itc_mean            = mean(I, 2);
    itc.itc_peak            = pk;
    itc.itc_peak_latency_ms = tPost(ip)' * 1000;
end
end
