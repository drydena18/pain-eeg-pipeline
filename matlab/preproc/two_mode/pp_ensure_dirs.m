function pp_ensure_dirs(cfg, subj)
%PP_ENSURE_DIRS pp_ensure_dirs(cfg) creates the experiment-level log folder;
%   pp_ensure_dirs(cfg, subj) also creates that participant's folders.
if nargin < 2
    P = pp_paths(cfg, 'sub-000');
    d = {P.dirs.logs};
else
    P = pp_paths(cfg, subj);
    d = struct2cell(P.dirs);
end
for i = 1:numel(d)
    if ~isfolder(d{i}), mkdir(d{i}); end
end
end