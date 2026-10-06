function m = spec_window_mask(timesSec, tmin, tmax, fs)
% SPEC_WINDOW_MASK  Inclusive [tmin, tmax] sample mask, quarter-sample tolerance
% V 1.0.0
%
% Python twin: src_spectral.src_window_mask. The tolerance stops a boundary
% sample being dropped by floating-point noise in the time axis (EEGLAB and
% MNE build it differently), so both pipelines select the same samples.
%
% Inputs:
%   timesSec : time axis in SECONDS (EEG.times / 1000)
%   tmin, tmax : window bounds in seconds
%   fs       : sampling rate in Hz
%
% Output:
%   m : logical row vector, same length as timesSec

tol = 0.25 / fs;
timesSec = double(timesSec(:))';
m = (timesSec >= tmin - tol) & (timesSec <= tmax + tol);
end
