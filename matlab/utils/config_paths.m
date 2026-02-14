function P = config_paths(exp_id, cfg)
% CONFIG_PATHS Centralized paths for CNED EEG pipeline (general across
% experiments)
% V 2.0.0
%
% Goals:
%   - Raw data can live on read-only mounts (/cifs/...)
%   - output/logs/QC live on writable location (local or writable share)
%   - exp_id selects experiment_specific folders via registry
%   - stable helpers for per-subject directories + naming
%
% Usage:
%   - P = config_paths('exp_01', cfg');
%   - P = config_paths(exp_id, cfg_from_jsonencode);
%
% Output:
%   P.RAW_ROOT
%   P.PROJ_ROOT
%   P.INPUT.EXP
%   P.RUN_ROOT
%   P.RESOURCE
%   P.CORE.ELP_FILE
%   P.SUBJECT_DIR(subjid)
%   P.STAGE_DIR(subjid, stageKey)   stageKey in {"EEG", "LOGS", "GC"}
%   P.NAMING.fname(subjid, tags, prefix)
%
% Notes:
%   - If you want a seperate folder per experiment within P.PROJ_ROOT, that
%   is exactly what P.RUN_ROOT does:
%       P.RUN_ROOT = fullfile(P.PROJ_ROOT, E.out_dirname, 'preproc')

% ----------------
% Args / defaults
% ----------------
if nargin < 1 || isempty(exp_id)
    exp_id = "exp01";
end
exp_id = string(exp_id);

if nargin < 2
    cfg = struct();
end

exp_id = string(exp_id);
exp_id = erase(exp_id, "_");

% --------------------
% Experiment Registry
% --------------------
R = experiment_registry();

if ~isfield(R, exp_id)
    error('config_paths:UnknownEXP', ...
        'Unknown exp_id "%s". Add it to experiment_registry() in config_paths.m', exp_id);
end

E = R.(exp_id);

% ----------------------
% Roots: RAW vs OUTPUTS
% ----------------------
% RAW root can be read-only
rawRootDefault = fullfile('/cifs', 'seminowicz', 'eegPainDatasets', 'CNED');

% Outputs must be writeable. Default to user home.
outRootDefault = fullfile(string(getenv("HOME")), "CNED_outputs");

% Allow JSON override:
% cfg.paths.raw_root, cfg.paths.proj_root or (out_root)
P = struct();

P.RAW_ROOT = get_cfg_string(cfg, ['paths', 'raw_root'], rawRootDefault);
P.PROJ_ROOT = get_cfg_string(cfg, ['paths', 'proj_root'], outRootDefault);
% Also support cfg.paths.out_root as alias
P.PROJ_ROOT = get_cfg_string(cfg, ['paths', 'out_root'], P.PROJ_ROOT);

% -------------------
% Raw + output paths
% -------------------
P.EXP = E; % keep registry entry for logging/debug

P.INPUT = struct();
P.INPUT.EXP = fullfile(P.RAW_ROOT, E.raw_dirname);

% All outputs for this experiment go under PROJ_ROOT/<out_dirname>/preproc
P.RUN_ROOT = fullfile(P.PROJ_ROOT, E.out_dirname, 'preproc');
P.RESOURCE = fullfile(P.RUN_ROOT, 'resources');

% ---------------
% Core Resources
% ---------------
P.CORE = struct();

% Default ELP lives in resources unless cfg overrides it
elpDefault = fullfile(P.RESOURCE, 'standard-10-5-cap385.elp');
P.CORE.ELP_FILE = elpDefault;

% If cfg.exp.channel_locs.use_elp and cfg.exp.channel_locs.elp_path_key
% exist, downstream can resolve the key via P.CORE.(...) or via your own
% key resovler

P.CORE.PARTICIPANTS_TSV = fullfile(P.RESOURCE, 'participants.tsv');
P.CORE.CSV_SINGLETRIAL = fullfile(P.RESOURCE, sprintf('participants_singletrial_%s.csv', E.raw_dirname));
P.CORE.DATASET_EVENTS_JSON = fullfile(P.INPUT.EXP, sprintf('task-%s_events.json', E.raw_dirname));

% --------------------------
% Per-Subject Folder Scheme
% --------------------------
P.SUBJECT_DIR = @(subjid) fullfile(P.RUN_ROOT, sprintf('sub-%03d', subjid));

P.STAGE = struct();
P.STAGE.FILTER      = '01_filter';
P.STAGE.NOTCH       = '02_notch';
P.STAGE.RESAMPLE    = '03_resample';
P.STAGE.REREF       = '04_reref';
P.STAGE.INITREJ     = '05_initrej';
P.STAGE.ICA         = '06_ica';
P.STAGE.EPOCH       = '07_epoch';
P.STAGE.BASE        = '08_base';
P.STAGE.LOGS        = 'LOGS';
P.STAGE.QC          = 'QC';

P.STAGE_DIR = @(subjid, stageKey) fullfile(P.SUBJECT_DIR(subjid), P.STAGE.(upper(string(stageKey))));

% --------------
% Naming Helper
% --------------
% Prefer cfg.exp.out_prefix if present; else fallback to registry; else
% fallback default
prefixDefault = 'analyzed_';
prefixFromCfg = get_cfg_string(cfg, ['exp', 'out_prefix'], '');
if strlength(prefixFromCfg) > 0
    defaultPrefix = prefixFromCfg;
elseif isfield(E, 'out_prefix') && strlength(string(E.out_prefix)) > 0
    defaultPrefix = string(E.out_prefix);
else
    defaultPrefix = prefixDefault;
end

P.NAMING = struct();
P.NAMING.default_prefix = defaultPrefix;
P.NAMING.fname = @(subjid, tags, prefix) local_fname(subjid, tags, prefix, defaultPrefix);

% -------------------------
% Ensure Output Dirs Exist
% -------------------------
% IMPORTANT: only mkdir under PROJ_ROOT (writable). Never touch RAW_ROOT.
ensure_dir(P.PROJ_ROOT);
ensure_dir(P.RUN_ROOT);
ensure_dir(P.RESOURCE);

end