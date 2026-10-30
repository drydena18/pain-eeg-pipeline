function y = pp_sanitize(x)
%PP_SANITIZE Recursively replace function handles with their text so a
%   config can be JSON-encoded and hashed.
if isa(x, 'function_handle')
    y = func2str(x);
elseif isstruct(x)
    y = x;
    f = fieldnames(x);
    for i = 1:numel(x)
        for j = 1:numel(f)
            y(i).(f{j}) = pp_sanitize(x(i).(f{j}));
        end
    end
elseif iscell(x)
    y = cellfun(@pp_sanitize, x, 'UniformOutput', false);
else
    y = x;
end
end