function val = get_cfg_string(cfg, keyPath, defaultVal)
% Safely get nested string from cfg struct:
% keyPath = ['paths', 'raw_root'] etc.
val = string(defaultVal);

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