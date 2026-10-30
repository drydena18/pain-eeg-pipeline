function h = pp_hash(x)
%PP_HASH SHA-256 hex digest of a numeric array, string, or struct.
%   Structs/cells are hashed via their JSON encoding, with function handles
%   converted to their text (pp_sanitize).
if isnumeric(x) || islogical(x)
    bytes = typecast(double(x(:)), 'uint8');
else
    bytes = unicode2native(jsonencode(pp_sanitize(x)), 'UTF-8');
end
md = java.security.MessageDigest.getInstance('SHA-256');
if ~isempty(bytes)
    md.update(typecast(bytes(:), 'int8'));
end
d = typecast(md.digest(), 'uint8');
h = sprintf('%02x', d);
end