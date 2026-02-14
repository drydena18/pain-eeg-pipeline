function s = defaultField(s, field, val)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = val;
end
end