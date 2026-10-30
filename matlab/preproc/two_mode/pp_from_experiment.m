function pp = pp_from_experiment(exp_id, subjects_override)
%PP_FROM_EXPERIMENT Two-mode runner config for ONE experiment.
%   Loads config/<exp_id>.json and prepares it exactly as expXX_preproc ->
%   preproc_default do (config_paths, the same defaults, the same subject
%   resolution), then wires the stage hooks from pp_hooks.m.
%
%   pp = pp_from_experiment('exp06');            % subjects in exp06.json
%   pp = pp_from_experiment('exp06', [1 4 17]);  % subset
%
%   Usually called through preproc_2mode('exp06', ...).

if nargin < 2, subjects_override = []; end
exp_id = char(erase(string(exp_id), "_"));

here = fileparts(mfilename('fullpath'));              % <repo>/matlab/preproc/two_mode
repo = fileparts(fileparts(fileparts(here)));
cfgPath = fullfile(repo, 'config', [exp_id '.json']);
if ~isfile(cfgPath)
    error('pp:exp', 'Config JSON not found: %s', cfgPath);
end
ecfg = load_cfg(cfgPath);
P    = config_paths(exp_id, ecfg);

% ---- Same validation / defaults as preproc_default (minus the dispatch) ----
mustHave(ecfg, 'exp', 'Missing cfg.exp in JSON.');
mustHave(ecfg, 'preproc', 'Missing cfg.preproc in JSON');
mustHave(ecfg.exp, 'out_prefix', 'Missing cfg.exp.out_prefix.');
if ~strcmp(char(string(ecfg.exp.id)), exp_id)
    error('pp:exp', '%s says exp.id = "%s", expected "%s".', cfgPath, string(ecfg.exp.id), exp_id);
end

if ~isempty(subjects_override)
    ecfg.exp.subjects = normalize_subject_ids(subjects_override);
elseif ~isfield(ecfg.exp, 'subjects') || isempty(ecfg.exp.subjects)
    tsv = {fullfile(string(P.INPUT.EXP), 'participants.tsv'), string(P.CORE.PARTICIPANTS_TSV)};
    tsv = tsv(cellfun(@(f) exist(f, 'file') == 2, tsv));
    if isempty(tsv)
        error('pp:exp', 'exp.subjects empty and participants.tsv not found in raw or resources.');
    end
    T = readtable(tsv{1}, 'FileType', 'text', 'Delimiter', '\t');
    ecfg.exp.subjects = normalize_subject_ids(extract_subject_column(T));
else
    ecfg.exp.subjects = normalize_subject_ids(ecfg.exp.subjects);
end
if isempty(ecfg.exp.subjects), error('pp:exp', 'No subjects resolved.'); end

ecfg.exp = defaultStruct(ecfg.exp, 'raw');
ecfg.exp.raw = defaultField(ecfg.exp.raw, 'pattern', '');
ecfg.exp.raw = defaultField(ecfg.exp.raw, 'search_recursive', true);
ecfg.exp = defaultStruct(ecfg.exp, 'channel_locs');
ecfg.exp.channel_locs = defaultField(ecfg.exp.channel_locs, 'use_elp', false);
ecfg.exp.channel_locs = defaultField(ecfg.exp.channel_locs, 'elp_path_key', "");
ecfg.exp = defaultStruct(ecfg.exp, 'montage');
ecfg.exp.montage = defaultField(ecfg.exp.montage, 'enabled', false);
ecfg.exp.montage = defaultField(ecfg.exp.montage, 'csv', "");
ecfg.exp.montage = defaultField(ecfg.exp.montage, 'select_ab_only', true);
ecfg.exp.montage = defaultField(ecfg.exp.montage, 'do_lookup', true);
ecfg.preproc = normalize_preproc_defaults(ecfg.preproc);

% ---- Runner config ----------------------------------------------------------
pp = pp_defaults();
pp.exp_id = exp_id;
pp.exp    = ecfg;      % stage hooks read the experiment config from cfg.exp
pp.P      = P;         % ... and config_paths output from cfg.P
pp.subjects     = arrayfun(@(k) sprintf('sub-%03d', k), ecfg.exp.subjects(:)', 'UniformOutput', false);
pp.paths.out    = char(P.RUN_ROOT);                     % <da_root>/<out_dirname>/preproc
pp.paths.prefix = char(string(ecfg.exp.out_prefix));
pp.paths.code   = repo;

pre = ecfg.preproc;
pp.mode1.avg_ref = logical(pre.reref.enabled) && strcmpi(string(pre.reref.mode), "average");
pp.mode1.reref_before_interp = true;                    % preproc_core: reref (04) before initrej (05)

% Copy the JSON blocks each mode depends on, so editing the JSON marks the
% right mode as out of date (scroll/prompt settings excluded: not used unattended).
m1 = rmfield(pre, intersect({'epoch', 'baseline', 'hilbert'}, fieldnames(pre)));
if isfield(m1, 'initrej'), m1.initrej = local_drop(m1.initrej, {'scroll'}); end
if isfield(m1, 'ica'), m1.ica = local_drop(m1.ica, {'iclabel', 'grid', 'scroll', 'confirm_scroll', 'train_scroll'}); end
pp.mode1.exp_preproc = m1;
pp.mode1.exp_input   = local_keep(ecfg.exp, {'raw', 'montage', 'channel_locs'});
if isfield(pre, 'ica')
    pp.mode2a.exp_review = local_keep(pre.ica, {'iclabel', 'grid'});
    if isfield(pre.ica, 'grid') && isfield(pre.ica.grid, 'n_ics'), pp.mode2a.n_show = pre.ica.grid.n_ics; end
end
pp.mode2b.exp_preproc = local_keep(pre, {'epoch', 'baseline', 'hilbert'});
if isfield(ecfg, 'features') && isfield(ecfg.features, 'bands') && isfield(ecfg.features.bands, 'slow_hz')
    pp.mode2b.exp_slow_hz = ecfg.features.bands.slow_hz;   % read by stage 09
end

% ---- Stage hooks (shared across experiments) ------------------------------
pp = pp_hooks(pp);

fprintf('[%s] two-mode: %d subject(s) | interp_policy = %s | badseg_policy = %s\n  -> %s\n', ...
    exp_id, numel(pp.subjects), pp.mode1.interp_policy, pp.mode1.badseg_policy, pp.paths.out);
end

% -------------------------------------------------------------------------
function s = local_keep(x, keys)
s = struct();
for i = 1:numel(keys)
    if isfield(x, keys{i}), s.(keys{i}) = x.(keys{i}); end
end
end

function x = local_drop(x, keys)
x = rmfield(x, intersect(keys, fieldnames(x)));
end