function featGA = spec_ga_mean_feat(featChan)
% SPEC_GA_MEAN_FEAT Channel-mean (omitnan) of a [nChan x nTr] feature struct
% V 1.0.0
%
% Reduces every field of a [nChan x nTr] feature struct to [nTr x 1] via
% mean across channels (dim 1, omitnan), returned as a column vector.
%
% Used for the pre_/post_/delta_ GA families in spectral_core.m.
%
% Input:
%   featChan : struct of [nChan x nTr] numeric fields
%
% Output:
%   featGA : struct of [nTr x 1] numeric fields, same field names

fn = fieldnames(featChan);
featGA = struct();
for i = 1:numel(fn)
    m = mean(featChan.(fn{i}), 1, 'omitnan');
    featGA.(fn{i}) = m(:);
end
end