function s = defaultStruct(s, field)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = struct();
end
end