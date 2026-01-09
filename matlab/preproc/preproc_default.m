function cfg = preproc_default(exp_id, P, cfg_in)
% PREPROC_DEFAULT Validate/normalize config then dispatch to preproc_core
% V 1.0.2
%
% Inputs:
%   exp_id  : string/char experiment id (e.g., "exp01")
%   P       : paths struct from config_paths(exp_id)
%   cfg_in  : struct from jsondecode(fileread(...))
%
% Output:
%   cfg     : normalized config (used by preproc_core)

% ---------------------
% Normalize exp_id type
% ---------------------
if ischar(exp_id) || isstring(exp_id)
    exp_id_str = string(exp_id);
else
    exp_id_str = exp_id;
end

cfg = cfg_in;

% ------------------------
% Basic required structure
% ------------------------
mustHave(cfg, 'exp', 'Missing cfg.exp in JSON.');
mustHave(cfg, 'preproc', 'Missing cfg.preproc in JSON.');

% Fill exp.id if absent
if ~isfield(cfg.exp, 'id') || isempty(cfg.exp.id)
    cfg.exp.id = exp_id_str; 
end

% Warn if mismatch (warn only; don't hard fail)
if isfield(cfg.exp, 'id') && ~isempty(cfg.exp.id)
    if string(cfg.exp.id) ~= exp_id_str
        warning('preproc_default:ExpIdMismatch', ...
            'exp_id=%s but cfg.exp.id=%s (continuing).', exp_id_str, string(cfg.exp.id));
    end
end

% Required naming prefix
mustHave(cfg.exp, 'out_prefix', 'Missing cfg.exp.out_prefix (e.g., "26BB_64_").');

% Ensure config_paths created naming struct
if ~isfield(P, 'NAMING') || ~isfield(P.NAMING, 'PREFIX')
    error('preproc_default:MissingNaming', 'P.NAMING.PREFIX missing in config_paths.m');
end

% Force the prefix to match JSON (char for sprintf stability)
P.NAMING.PREFIX = char(string(cfg.exp.out_prefix));

% -------------------------------------------
% Subjects: JSON override or participants.tsv
% -------------------------------------------
if ~isfield(cfg.exp, 'subjects') || isempty(cfg.exp.subjects)

    % Resolve participants.tsv path (support both nested and flat styles)
    tsvPath = '';
    if isfield(P, 'CORE')
        if isfield(P.CORE, 'PARTICIPANTS') && isfield(P.CORE.PARTICIPANTS, 'TSV')
            tsvPath = P.CORE.PARTICIPANTS.TSV;
        elseif isfield(P.CORE, 'PARTICIPANTS_TSV')
            tsvPath = P.CORE.PARTICIPANTS_TSV;
        end
    end

    if isempty(tsvPath) || ~exist(tsvPath, 'file')
        error('preproc_default:MissingParticipants', ...
            'cfg.exp.subjects empty and participants.tsv not found (check P.CORE.*): %s', string(tsvPath));
    end

    % Read participants table
    T = readtable(tsvPath, 'FileType', 'text', 'Delimiter', '\t');

    % Try common column names
    vn = lower(string(T.Properties.VariableNames));

    if any(vn == "subjid")
        subs_raw = T.(T.Properties.VariableNames{find(vn == "subjid", 1)});
    elseif any(vn == "participant_id")
        subs_raw = T.(T.Properties.VariableNames{find(vn == "participant_id", 1)});
    elseif any(vn == "subject")
        subs_raw = T.(T.Properties.VariableNames{find(vn == "subject", 1)}); 
    else
        error('preproc_default:BadParticipantsTSV', ...
            'participants.tsv needs a column like subjid / participant_id / subject.');
    end

    cfg.exp.subjects = normalize_subject_ids(subs_raw);
else
    cfg.exp.subjects = normalize_subject_ids(cfg.exp.subjects);
end

if isempty(cfg.exp.subjects)
    error('preproc_default:NoSubjects', ...
        'No subjects resolved. Check cfg.exp.subjects or participants.tsv.');
end

% -----------------------------
% Raw input settings (defaults)
% -----------------------------
cfg.exp.raw = defaultStruct(cfg.exp, 'raw');
cfg.exp.raw = defaultField(cfg.exp.raw, 'pattern', '');
cfg.exp.raw = defaultField(cfg.exp.raw, 'search_recursive', true);

% Channel location settings (defaults)
cfg.exp.channel_locs = defaultStruct(cfg.exp, 'channel_locs');
cfg.exp.channel_locs = defaultField(cfg.exp.channel_locs, 'use_elp', false);

% Montage defaults (handled in preproc_core after import)
cfg.exp = defaultStruct(cfg.exp, 'montage');
cfg.exp.montage = defaultField(cfg.exp.montage, 'enabled', false);
cfg.exp.montage = defaultField(cfg.exp.montage, 'csv', "");
cfg.exp.montage = defaultField(cfg.exp.montage, 'select_ab_only', true);
cfg.exp.montage = defaultField(cfg.exp.montage, 'do_lookup', true);

% ----------------------------------
% Preproc defaults + minimal validation
% ----------------------------------
cfg.preproc = normalize_preproc_defaults(cfg.preproc);

% --------
% Dispatch
% --------
preproc_core(P, cfg); % subjects_override handled inside preproc_core if you pass it there

end

% -------------------
% Helpers
% -------------------
function mustHave(s, field, msg)
if ~isfield(s, field)
    error('preproc_default:MissingField', '%s', msg);
end
end

function s = defaultStruct(s, field)
% Ensure s.(field) exists and is a struct
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = struct();
end
end

function s = defaultField(s, field, val)
% Ensure s.(field) exists and is non-empty
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = val;
end
end

function subs = normalize_subject_ids(subs_raw)
% Convert subjects into a numeric column vector.
% Accepts:
%   [1 2 3]
%   {"001","002"}
%   ["sub-001","sub-002"]
%   ["001","002"]
% Returns:
%   [1;2;3]

if isnumeric(subs_raw)
    subs = subs_raw(:);
    return;
end

subs_str = string(subs_raw);
subs = nan(numel(subs_str), 1);

for i = 1:numel(subs_str)
    x = subs_str(i);
    % Extract trailing digits (e.g., "sub-001" -> "001")
    tok = regexp(char(x), '(\d+)$', 'tokens', 'once'); 
    if ~isempty(tok)
        subs(i) = str2double(tok{1});
    end
end

if any(isnan(subs))
    badIdx = find(isnan(subs), 1, 'first');
    error('preproc_default:BadSubjectId', ...
        'Could not parse subject id from value "%s".', subs_str(badIdx));
end

subs = subs(:);
end

function Pp = normalize_preproc_defaults(Pp)
% Ensures each preproc block has enabled/tag + required params with defaults

% FILTER
Pp = ensureBlock(Pp, 'filter', true, 'fir');
Pp.filter = defaultField(Pp.filter, 'type', 'fir');
Pp.filter = defaultField(Pp.filter, 'highpass_hz', 0.5);
Pp.filter = defaultField(Pp.filter, 'lowpass_hz', 40);

% NOTCH
Pp = ensureBlock(Pp, 'notch', true, 'notch60');
Pp.notch = defaultField(Pp.notch, 'freq_hz', 60);
Pp.notch = defaultField(Pp.notch, 'bw_hz', 2);

% RESAMPLE
Pp = ensureBlock(Pp, 'resample', false, 'rs500');
Pp.resample = defaultField(Pp.resample, 'target_hz', []);

% REREF
Pp = ensureBlock(Pp, 'reref', true, 'reref');
Pp.reref = defaultField(Pp.reref, 'mode', 'average');
Pp.reref = defaultField(Pp.reref, 'channels', []);

% INITREJ
Pp = ensureBlock(Pp, 'initrej', false, 'initrej');
Pp.initrej = defaultStruct(Pp.initrej, 'badchan');
Pp.initrej = defaultStruct(Pp.initrej, 'badseg');
Pp.initrej.badchan = defaultField(Pp.initrej.badchan, 'enabled', false);
Pp.initrej.badseg  = defaultField(Pp.initrej.badseg,  'enabled', false);

% ICA
Pp = ensureBlock(Pp, 'ica', true, 'ica');
Pp.ica = defaultField(Pp.ica, 'method', 'runica');

% ICLabel sub-block
Pp.ica = defaultStruct(Pp.ica, 'iclabel');
Pp.ica.iclabel = defaultField(Pp.ica.iclabel, 'enabled', false);
Pp.ica.iclabel = defaultField(Pp.ica.iclabel, 'auto_reject', false);
Pp.ica.iclabel = defaultStruct(Pp.ica.iclabel, 'thresholds');

% EPOCH
Pp = ensureBlock(Pp, 'epoch', true, 'epoch');
Pp.epoch = defaultField(Pp.epoch, 'event_types', {});
Pp.epoch = defaultField(Pp.epoch, 'tmin_sec', -1.0);
Pp.epoch = defaultField(Pp.epoch, 'tmax_sec', 2.0);

% BASELINE
Pp = ensureBlock(Pp, 'baseline', true, 'base');
Pp.baseline = defaultField(Pp.baseline, 'window_sec', [-0.5 0]);

% Minimal validation
if Pp.epoch.enabled && isempty(Pp.epoch.event_types)
    error('preproc_default:EpochMissingEvents', ...
        'cfg.preproc.epoch.enabled is true, but epoch.event_types is empty.');
end
end

function S = ensureBlock(S, field, defaultEnabled, defaultTag)
if ~isfield(S, field) || isempty(S.(field))
    S.(field) = struct();
end
S.(field) = defaultField(S.(field), 'enabled', defaultEnabled);
S.(field) = defaultField(S.(field), 'tag', defaultTag);
end