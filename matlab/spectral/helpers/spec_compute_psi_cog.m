function feat = spec_compute_psi_cog(feat)
% SPEC_COMPUTE_PSI_COG Add the BI x (CoG - 10) interaction term to a feature struct
% V 1.0.0
%
% psu_cog = sf_balance .* (paf_cog_hz - 10)
%
% Inputs:
%   feat : struct, must contain .sf_balance and .paf_cog_hz
%          (both [nChan x nTr] or [nTr x 1] / [1 x nTr] for a GA struct)
%
% Output:
%   feat : same struct, with .psi_cog added

if ~isfield(feat, 'sf_balance') || ~isfield(feat, 'paf_cog_hz')
    error('spec_compute_psi_cog:MissingFields', ...
        'feat must contain sf_balance and paf_cog_hz (got: %s).', ...
        strjoin(fieldnames(feat), ', '));
end

feat.psi_cog = feat.sf_balance .* (feat.paf_cog_hz - 10);
end