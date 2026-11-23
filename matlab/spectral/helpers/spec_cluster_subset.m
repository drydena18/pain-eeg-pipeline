function [members, Csub] = spec_cluster_subset(C, nChan)
% SPEC_CLUSTER_SUBSET  Electrodes needed by the clusters, and indices into them
% V 1.0.0
%
% Used by spectral_core when spectral.units = "clusters": only the member
% electrodes go through spec_compute_subject_features, so cluster indices
% must point into that subset instead of the full montage.
%
% Inputs:
%   C     : struct array from spec_resolve_clusters (idx into the full montage)
%   nChan : number of channels in the full montage
%
% Outputs:
%   members : sorted unique electrode indices used by ANY cluster
%             (including clusters below min_channels, so the QC and the
%             feature computation see the same electrodes)
%   Csub    : C with .idx re-expressed as indices into members

members = [];
for c = 1:numel(C)
    members = [members; C(c).idx(:)]; %#ok<AGROW>
end
members = unique(members);

map = zeros(nChan, 1);
map(members) = 1:numel(members);
Csub = C;
for c = 1:numel(C)
    Csub(c).idx = reshape(map(C(c).idx), size(C(c).idx));
end
end
