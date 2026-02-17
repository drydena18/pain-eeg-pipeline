function S = ensureBlock(S, field, defaultEnabled, defaultTag)
if ~isfield(S, field) || isempty(S.(field))
    S.(field) = struct();
end
S.(field) = defaultField(S.(field), 'enabled', defaultEnabled);
S.(field) = defaultField(S.(field), 'tag', defaultTag);
end