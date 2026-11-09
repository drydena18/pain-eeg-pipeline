function T = pp_manifest(cfg, subs, errs)
%PP_MANIFEST One row per participant: state, rank, interpolated channels,
%   ICs rejected, and any error from this run. Written to logs/pp_manifest.csv.
subs = cellstr(subs(:));
n = numel(subs);
if nargin < 3, errs = repmat({''}, n, 1); end

status  = cell(n, 1);
detail  = repmat({''}, n, 1);
stale   = false(n, 1);
rnk     = nan(n, 1);
nInterp = nan(n, 1);
nRej    = nan(n, 1);
for i = 1:n
    [status{i}, info] = pp_status(cfg, subs{i});
    stale(i) = info.ica_stale;
    detail{i} = info.reason;
    if ~isempty(info.ica)
        rnk(i)     = info.ica.rank;
        nInterp(i) = numel(info.ica.interp_labels);
    end
    if ~isempty(info.decisions)
        nRej(i) = numel(info.decisions.rejected);
    end
end

T = table(subs, status, stale, rnk, nInterp, nRej, detail, errs(:), 'VariableNames', ...
    {'subject', 'status', 'ica_stale', 'rank', 'n_interp', 'n_ic_rejected', 'detail', 'last_error'});
P = pp_paths(cfg, 'x');
writetable(T, P.manifest);
end
