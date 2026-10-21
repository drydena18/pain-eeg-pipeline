function val = get_cfg_string(cfg, keyPath, defaultVal)
% GET_CFG_STRING Safely get a nested string from a cfg struct.
% V 1.1.0
%
%   val = get_cfg_string(cfg, {'paths', 'da_root'}, defaultVal)
%
% keyPath must be a CELL array (or string array) of field names. A char
% concatenation such as ['paths', 'da_root'] is ONE char vector
% ('pathsda_root'), not two keys.
%
% V1.1.0: a char keyPath is treated as a single field name. Before, it was
% walked one CHARACTER at a time ('p', 'a', ...), so every lookup written as
% ['a', 'b'] silently returned the default (this is why config_paths.m never
% read cfg.paths / cfg.exp.out_prefix).

val = string(defaultVal);

if ischar(keyPath)
    keyPath = {keyPath};
end

try
    cur = cfg;
    for i = 1:numel(keyPath)
        k = char(keyPath(i));
        if ~isstruct(cur) || ~isfield(cur, k)
            return;
        end
        cur = cur.(k);
    end
    if isstring(cur) || ischar(cur)
        s = string(cur);
        if strlength(s) > 0
            val = s;
        end
    end
catch
    % fall back to defaultVal
end
end