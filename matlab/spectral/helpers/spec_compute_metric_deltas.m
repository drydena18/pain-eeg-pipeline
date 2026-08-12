function deltaFeat = spec_compute_metric_deltas(featPre, featPost)
% SPEC_COMPUTE_METRIC_DELTAS delta_<name> = post_<name> - pre_<name>
% V 1.1.0
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

if nargin < 2
    error('spec_compute_metric_deltas:MissingArg', ...
        'spec_compute_metric_deltas requires two arguments (featPre, featPost); got %d.', nargin);
end

if ~isstruct(featPre) || isempty(featPre)
    error('spec_compute_metric_deltas:BadPreFeat', ...
        'preFeat must be a non-empty struct; got class "%s", isempty = %d.', class(featPre), isempty(featPost));
end

if ~isstruct(featPost) || isempty(featPost)
    error('spec_compute_metric_deltas:BadPostFeat', ...
        'postFeat must be a non-empty struct; got class "%s", isempty = %d.', class(featPost), isempty(featPost));
end

fnPre = fieldnames(featPre);
fnPost = fieldnames(featPost);
common = intersect(fnPre, fnPost);
if isempty(common)
    error('spec_compute_metric_deltas:NoCommonFields', ...
        'preFeat and postFeat share no field names.\n preFeat fields: %s\n postFeat fields: %s', ...
        strjoin(fnPre, ', '), strjoin(fnPost, ', '));
end

onlyPre = setdiff(fnPre, fnPost);
onlyPost = setdiff(fnPost, fnPre);
if ~isempty(onlyPre)
    warning('spec_compute_metric_deltas:PreOnly', ...
        'Fields present in postFeat but not preFeat (skipped): %s', strjoin(onlyPre, ', '));
end
if ~isempty(onlyPost)
    warning('spec_compute_metric_deltas:PostOnly', ...
        'Fields present in preFeat but not postFeat (skipped): %s', strjoin(onlyPost, ', '));
end

deltaFeat = struct();
for i = 1:numel(common)
    fn = common{i};
    try
        deltaFeat.(['delta_', fn]) = featPost.(fn) - featPre.(fn);
    catch ME
        error('spec_compute_metric_deltas:SubtractFailed', ...
            'Failed computing delta_%s = featPost.%s - featPre.%s (sizes: post = %s, pre = %s). Underlying error: %s.', ...
            fn, fn, fn, mat2str(size(featPost.(fn))), mat2str(size(featPre.(fn))), ME.message);
    end
end
end