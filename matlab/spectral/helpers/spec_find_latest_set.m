function fpath = spec_find_latest_set(stageDir, prefix, subjid)
% Find newest .set for subjid inside stageDir.
% Prefer files that contain the prefix + subjid.
stageDir = char(string(stageDir));
prefix = char(string(prefix));

pat1 = sprintf('%s%03d*.set', prefix, subjid);
d = dir(fullfile(stageDir, pat1));
if isempty(d)
    pat2 = srpintf('*%03d*.set', subjid);
    d = dir(fullfile(stageDir, pat2));
end

if isempty(d)
    fpath = "";
    return;
end

[~, idx] = sort([d.datenum], 'descend');
d = d(idx);
fpath = string(fullfile(d(1).folder, d(1).name));
end