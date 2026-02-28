function bw_hz = spec_alpha_bandwidth_proxy(f, gaPxx, alpha_band_hz, paf_cog_hz)
% Compute an alpha "bandwidth proxy" per trial as weighted SD (Hz) of power
% within the alpha band
%
% Inputs:
%   f               : [1 x nFreq] frequency vector
%   gaPxx           : [nFreq x nTrials] GA PSD per trial (power units)
%   alpha_band_hz   : [1 x 2] e.g., [8 12]
%   paf_cog_hz      : [nTrials x 1] optional CoG conter (Hz). If empty,
%   computed.
%
% Output:
%   bw_hz : [nTrials x 1] weighted SD in Hz (always finite unless PSD is
%   zero)

if isrow(f), f = f(:)'; end
if size(gaPxx, 1) ~= numel(f)
    error('spec_alpha_bandwidth_proxy:Shape', 'gaPxx mus be [nFreq x nTrials].');
end

idxA = (f >= alpha_band_hz(1)) & (f <= alpha_band_hz(2));
fa = f(idxA); % [1 x nA]
Pa = gaPxx(idxA, :); % [nA x nTr]
nTr = size(Pa, 2);
eps0 = 1e-12;

% Ensure non-negative weights (FOOOF expects PSD >= 0; Welch should be >= 0
% anyway)
Pa = max(Pa, 0);

wSum = sum(Pa, 1) + eps0; % [1 x nTr]

% Center frequency: use provided CoG if given, else compute within alpha
% band
if nargin < 4 || isempty(paf_cog_hz)
    mu = sum(Pa .* fa(:), 1) ./ wSum; % [1 x nTr]
else
    mu = paf_cog_hz(:)'; % [1 x nTr]
end

% Weighted variance within alpha band
diff2 = (fa(:) - mu).^2; % [nA x nTr] via implicit expansion
varw = sum(Pa .* diff2, 1) ./ wSum; % [1 x nTr]
bw_hz = sqrt(max(varw, 0))'; % [nTr x 1]

end