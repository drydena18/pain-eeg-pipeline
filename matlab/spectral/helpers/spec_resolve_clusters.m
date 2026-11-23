function C = spec_resolve_clusters(ccfg, chanLabels)
% SPEC_RESOLVE_CLUSTERS  Turn cfg.spectral.clusters into channel indices
% V 1.0.0
%
% Two kinds of cluster definition (both optional, any number of each):
%
%   lateral.<Region>.left / .right : one electrode list per hemisphere.
%       Each region becomes TWO units, named relative to the stimulated
%       hand (ccfg.stim_side):
%           <Region>_contra  = hemisphere opposite the stimulated hand
%           <Region>_ipsi    = hemisphere on the same side
%       e.g. stim_side = "right" -> S1_contra = lateral.S1.left (C3 side).
%
%   fixed.<Name> : one electrode list, used as given (e.g. a midline or
%       bilateral cluster such as C3/Cz/C4).
%
% Matching is case-insensitive and ignores surrounding whitespace, so
% "CZ" (26ByBiosemi relabel) matches "Cz". Requested electrodes that are not
% in the recording are reported in .missing; a cluster with fewer than
% ccfg.min_channels electrodes found is marked ok = false (all-NaN rows).
%
% Inputs:
%   ccfg       : cfg.spectral.clusters (normalised by spectral_default)
%   chanLabels : channel labels of the loaded EEG (cellstr or string array)
%
% Output:
%   C : struct array, one element per output unit, fields
%       name      char      e.g. 'S1_contra'
%       kind      char      'lateral' or 'fixed'
%       requested cellstr   electrodes asked for, in config order
%       idx       double    indices into chanLabels of those found
%       found     cellstr   labels as they appear in the data
%       missing   cellstr   requested but not in the data
%       ok        logical   numel(idx) >= min_channels

labs  = local_cellstr(chanLabels);
labsU = upper(strtrim(labs));

side = lower(strtrim(char(ccfg.stim_side)));
switch side
    case 'right', contra = 'left';
    case 'left',  contra = 'right';
    otherwise
        error('spec_resolve_clusters:BadStimSide', ...
            'clusters.stim_side must be "left" or "right" (got "%s").', side);
end
minCh = double(ccfg.min_channels);

defs = cell(0, 3);   % {name, kind, list}

if isfield(ccfg, 'lateral') && isstruct(ccfg.lateral)
    regs = local_fields(ccfg.lateral);
    for r = 1:numel(regs)
        R = ccfg.lateral.(regs{r});
        if ~(isstruct(R) && isfield(R, 'left') && isfield(R, 'right'))
            error('spec_resolve_clusters:BadLateral', ...
                'clusters.lateral.%s needs both "left" and "right" electrode lists.', regs{r});
        end
        defs(end+1, :) = {[regs{r} '_contra'], 'lateral', R.(contra)}; %#ok<AGROW>
        defs(end+1, :) = {[regs{r} '_ipsi'],   'lateral', R.(side)};   %#ok<AGROW>
    end
end

if isfield(ccfg, 'fixed') && isstruct(ccfg.fixed)
    names = local_fields(ccfg.fixed);
    for k = 1:numel(names)
        defs(end+1, :) = {names{k}, 'fixed', ccfg.fixed.(names{k})}; %#ok<AGROW>
    end
end

nC = size(defs, 1);
C = repmat(struct('name', '', 'kind', '', 'requested', {{}}, 'idx', [], ...
                  'found', {{}}, 'missing', {{}}, 'ok', false), nC, 1);
if nC == 0
    return;
end

% Unit names must be unique and must not shadow a real channel label,
% because cluster rows are appended to the channel CSVs.
nameU = upper(defs(:, 1));
if numel(unique(nameU)) ~= numel(nameU)
    error('spec_resolve_clusters:DuplicateName', 'Cluster names are not unique: %s', ...
        strjoin(defs(:, 1)', ', '));
end
clash = intersect(nameU, labsU);
if ~isempty(clash)
    error('spec_resolve_clusters:NameClash', ...
        'Cluster name(s) %s equal a channel label; rename the cluster.', strjoin(clash(:)', ', '));
end

for c = 1:nC
    req  = local_cellstr(defs{c, 3});
    reqU = upper(strtrim(req));
    [tf, loc] = ismember(reqU, labsU);
    C(c).name      = defs{c, 1};
    C(c).kind      = defs{c, 2};
    C(c).requested = req;
    C(c).idx       = loc(tf);
    C(c).found     = labs(loc(tf));
    C(c).missing   = req(~tf);
    C(c).ok        = numel(C(c).idx) >= minCh;
end
end

% -------------------------------------------------------------------------
function out = local_cellstr(x)
% Row cellstr from char, cellstr, or string array (jsondecode gives a
% column cellstr for a JSON array of strings, or a char for one string).
if ischar(x)
    out = {x};
elseif iscell(x)
    out = cellfun(@(v) char(v), x(:)', 'UniformOutput', false);
else
    out = cellstr(x);   % string array (MATLAB)
    out = out(:)';
end
end

function fn = local_fields(s)
% Field names, skipping JSON comment keys ("_note" -> x_note after jsondecode).
fn = fieldnames(s);
keep = ~(strncmp(fn, 'x_', 2) | strncmp(fn, '_', 1));
fn = fn(keep);
end