function cfg = preproc_default(P, cfg_in, subjects_override)
% PREPROC_DEFAULT Validate/normalize config then dispatch to preproc_core
% V 1.1.0
%
% Inputs:
%   P                 : paths struct from config_paths(exp_id, cfg)
%   cfg_in            : struct from jsondecode(fileread(...))
%   subjects_override : [] to use cfg.exp.subjects; otherwise numeric
%   vector
%
% Output:
%   cfg     : normalized config (passed into preproc_core)

if nargin < 1 || isempty(P)
    error('preproc_default:MissingP', 'P is required.');
end
if nargin < 2 || isempty(cfg_in)
    error('preproc_default:MissingCfg', 'cfg_in is required.');
end
if nargin < 3
    subjects_override = [];
end

cfg = cfg_in;

% -------------------------
% Basic required structure 
% -------------------------
mustHave(cfg, 'exp', 'Missing cfg.exp in JSON.');
mustHave(cfg, 'preproc', 'Missing cfg.preproc in JSON');

% Fill exp.id from P if absent
if ~isfield(cfg.exp, 'id') || isempty(cfg.exp.id)
    if isfield(P, 'EXP') && isfield(P.EXP, 'id')
        cfg.exp.id = string(P.EXP.id);
    else
        cfg.exp.id = 'unknown_exp';
    end
end

% --------------------
% Naming / out_prefix
% --------------------
mustHave(cfg.exp, 'out_prefix', 'Missing cfg.exp.out_prefix (e.g., "26BB_62_").');

% Ensure config_paths naming exists
if ~isfield(P, 'NAMING') || ~isfield(P.NAMING, 'fname')
    error('preproc_default:MissingNaming', 'P.NAMING.fname missing in config_paths.m');
end

% Override naming default prefix from JSON (for consistent filenames)
P.NAMING.default_prefix = string(cfg.exp.out_prefix);

% --------------------------------
% Subjects: override > JSON > TSV
% --------------------------------
if ~isempty(subjects_override)
    cfg.exp.subjects = normalize_subject_ids(subjects_override);
elseif ~isfield(cfg.exp, 'subjects') || isempty(cfg.exp.subjects)
    % Try raw participants.tsv first, then resources
    tsvCandidates = {};

    % raw dataset participants.tsv (common BIDS-ish)
    if isfield(P, 'INPUT') && isfield(P.INPUT, 'EXP')
        tsvCandidates{end+1} = fullfile(string(P.INPUT.EXP), 'participants.tsv');
    end

    % resources participants.tsv (config_paths sets this)
    if isfield(P, 'CORE') && isfield(P.CORE, 'PARTICIPANTS_TSV')
        tsvCandidates{end+1} = string(P.CORE.PARTICIPANTS_TSV);
    end

    tsvPath = '';
    for i = 1:numel(tsvCandidates)
        if exist(tsvCandidates{i}, 'file')
            tsvPath = string(tsvCandidates{i});
            break;
        end
    end

    if strlength(tsvPath) == 0
        error('preproc_default:MissingParticipants', ...
            'cfg.exp.subjects empty and participants.tsv not found in raw or resources.');
    end

    T = readtable(tsvPath, 'FileType', 'text', 'Delimiter', '\t');
    cfg.exp.subjects = normalize_subject_ids(extract_subject_column(T));
else
    cfg.exp.subjects = normalize_subject_ids(cfg.exp.subjects);
end

if isempty(cfg.exp.subjects)
    error('preproc_default:NoSubjects', 'No subjects resolved.');
end

% ------------------------------
% Raw Input Settings (defaults)
% ------------------------------
cfg.exp = defaultStruct(cfg.exp, 'raw');
cfg.exp.raw = defaultField(cfg.exp.raw, 'pattern', '');
cfg.exp.raw = defaultField(cfg.exp.raw, 'search_recursive', true);

% Channel location settings (defaults)
cfg.exp = defaultStruct(cfg.exp, 'channel_locs');
cfg.exp.channel_locs = defaultField(cfg.exp.channel_locs, 'use_elp', false);
cfg.exp.channel_locs = defaultField(cfg.exp.channel_locs, 'elp_path_key', "");

% Montage defaults
cfg.exp = defaultStruct(cfg.exp, 'montage');
cfg.exp.montage = defaultField(cfg.exp.montage, 'enabled', false);
cfg.exp.montage = defaultField(cfg.exp.montage, 'csv', "");
cfg.exp.montage = defaultField(cfg.exp.montage, 'select_ab_only', true);
cfg.exp.montage = defaultField(cfg.exp.montage, 'do_lookup', true);

% --------------------------------------
% Preproc defaults + minimal validation
% --------------------------------------
cfg.preproc = normalize_preproc_defaults(cfg.preproc);

% ---------------
% Log run header
% ---------------
fprintf('[%s] Subjects (%d): %s\n', string(cfg.exp.id), numel(cfg.exp.subjects), mat2str(cfg.exp.subjects(:)'));

% ----------
% Dispatch
% ----------
preproc_core(P, cfg);

end