function ga = spec_squeeze_ga_features(featGA)
% featGA fields are [1 x trial] or [trial x 1]; normalize to [trial x 1]
ga = struct();
fn = fieldnames(featGA);
for i = 1:numel(fn)
    x = featGA.(fn{i});
    x = squeeze(x);
    if isrow(x), x = x'; end
    ga.(fn{i}) = x;
end
end