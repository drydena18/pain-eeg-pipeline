function deltaFeat = spec_compute_metric_deltas(featPre, featPost)
% SPEC_COMPUTE_METRIC_DELTAS delta_<name> = post_<name> - pre_<name>
% V 1.0.0
%
% Generic raw pre->post-change, computed for every field common to both
% UNPREFIXED preFeat structs. This is a plain different for every metric,
% regardless of whether the metric is a bounded ratio (sf_balance, psi_cog, ...)
% or a power quantity (pow_slow_alpha, ...).
%
% Inputs:
%   preFeat, postFeat : structs with matching, unprefixed field names,
%                       each field [nChan x nTr] (or [nTr x 1] for GA)
%
% Output:
%   deltaFeat : strict with fields delta_<name>

fnPre = fieldnames(featPre);
fnPost = fieldnames(featPost);
common = intersect(fnPre, fnPost);

onlyPre = setdiff(fnPre, fnPost);
onlyPost = setdiff(fnPost, fnPre);
if ~isempty(onlyPre)
    warning('spec_compute_metric_deltas:PreOnly', ...
        'Fields present in postFeat but not preFeat (skipped): %s', strjoin(onlyPost, ', '));
end

deltaFeat = struct();
for i = 1:numel(common)
    fn = common{i};
    deltaFeat.(['delta_', fn]) = postFeat.(fn) - preFeat.(fn);
end
end