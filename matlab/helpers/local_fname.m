function fname = local_fname(subjid, tags, prefix, defaultPrefix)
% fname = <prefix><subjid>_<tag>_<tag>.set
if nargin < 3 || isempty(prefix)
    prefix = defaultPrefix;
end
prefix = string(prefix);

if nargin < 2 || isempty(tags)
    tagStr = "";
else
    if isstring(tags), tags = cellstr(tags); end
    if ischar(tags), tags = {tags}; end
    tags = cellfun(@(x) string(x), tags, 'UniformOutput', false);
    tagStr = "_" + strjoin([tags{:}], "_");
end

fname = sprintf('%s%03d%s.set', prefix, subjid, tagStr);
fname = char(fname);
end
