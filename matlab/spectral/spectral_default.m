function cfg = spectral_default(P, cfg_in, subjects_override)
% SPECTRAL_DEFAULT Validate/normalise cfg.spectral then dispatch to spectral_core
% V 1.1.0
%
% New in V1.1:
%   cfg.spectral.windows   -- pre/post-stim windows for interaction metrics
%   cfg.spectral.lep       -- LEP extraction config
%   cfg.spectral.phase     -- Hilbert phase config
%   cfg.spectral.trial_spectral  -- now fully defaulted (was partial)
%
% Inputs:
%   P                 : paths struct from config_paths(exp_id, cfg)
%   cfg_in            : struct from load_cfg(json)
%   subjects_override : [] to use cfg.exp.subjects; otherwise numeric vector

if nargin < 1 || isempty(P)
    error('spectral_default:MissingP', 'P is required.');
end
if nargin < 2 || isempty(cfg_in)
    error('spectral_default:MissingCfg', 'cfg_in is required.');
end
if nargin < 3
    subjects_override = [];
end

cfg = cfg_in;

mustHave(cfg, 'exp',      'Missing cfg.exp in JSON.');
mustHave(cfg, 'spectral', 'Missing cfg.spectral in JSON.');

if ~isfield(cfg.exp, 'id') || isempty(cfg.exp.id)
    if isfield(P, 'EXP') && isfield(P.EXP, 'id')
        cfg.exp.id = string(P.EXP.id);
    else
        cfg.exp.id = "unknown_exp";
    end
else
    cfg.exp.id = string(cfg.exp.id);
end

mustHave(cfg.exp, 'out_prefix', 'Missing cfg.exp.out_prefix (e.g., "26BB_62_").');

% ------------------
% Resolve subjects
% ------------------
if ~isempty(subjects_override)
    cfg.exp.subjects = normalize_subject_ids(subjects_override);

elseif ~isfield(cfg.exp, 'subjects') || isempty(cfg.exp.subjects)
    tsvCandidates = {};
    if isfield(P, 'INPUT') && isfield(P.INPUT, 'EXP')
        tsvCandidates{end+1} = fullfile(string(P.INPUT.EXP), 'participants.tsv');
    end
    if isfield(P, 'CORE') && isfield(P.CORE, 'PARTICIPANTS_TSV')
        tsvCandidates{end+1} = string(P.CORE.PARTICIPANTS_TSV);
    end

    tsvPath = "";
    for i = 1:numel(tsvCandidates)
        if exist(tsvCandidates{i}, 'file')
            tsvPath = string(tsvCandidates{i});
            break;
        end
    end

    if strlength(tsvPath) == 0
        error('spectral_default:MissingParticipants', ...
            'cfg.exp.subjects empty and participants.tsv not found in raw or resources.');
    end

    T = readtable(tsvPath, 'FileType', 'text', 'Delimiter', '\t');
    cfg.exp.subjects = normalize_subject_ids(extract_subject_column(T));
else
    cfg.exp.subjects = normalize_subject_ids(cfg.exp.subjects);
end

if isempty(cfg.exp.subjects)
    error('spectral_default:NoSubjects', 'No subjects resolved.');
end

% ---------------------------
% Normalise spectral config
% ---------------------------
cfg.spectral = defaultField(cfg.spectral, 'enabled',     true);
cfg.spectral = defaultField(cfg.spectral, 'input_stage', "08_base");

% ---- PSD params ----
cfg.spectral = defaultStruct(cfg.spectral, 'psd');
cfg.spectral.psd = defaultField(cfg.spectral.psd, 'fmin_hz',      1);
cfg.spectral.psd = defaultField(cfg.spectral.psd, 'fmax_hz',      80);
cfg.spectral.psd = defaultField(cfg.spectral.psd, 'window_sec',   2.0);
cfg.spectral.psd = defaultField(cfg.spectral.psd, 'overlap_frac', 0.5);
cfg.spectral.psd = defaultField(cfg.spectral.psd, 'nfft',         0);

% ---- Alpha band ----
cfg.spectral = defaultStruct(cfg.spectral, 'alpha');
cfg.spectral.alpha = defaultField(cfg.spectral.alpha, 'alpha_hz', [8 12]);
cfg.spectral.alpha = defaultField(cfg.spectral.alpha, 'slow_hz',  [8 10]);
cfg.spectral.alpha = defaultField(cfg.spectral.alpha, 'fast_hz',  [10 12]);

% ---- Pre/post-stim windows (NEW) ----
% pre_sec  : baseline for BI_pre, LR_pre, CoG_pre, and ERD pre-power
% post_sec : response for ERD post-power
% Both in seconds relative to stimulus onset (t = 0).
cfg.spectral = defaultStruct(cfg.spectral, 'windows');
cfg.spectral.windows = defaultField(cfg.spectral.windows, 'pre_sec',  [-1.0, -0.1]);
cfg.spectral.windows = defaultField(cfg.spectral.windows, 'post_sec', [ 0.1,  0.8]);

% ---- LEP extraction (NEW) ----
cfg.spectral = defaultStruct(cfg.spectral, 'lep');
cfg.spectral.lep = defaultField(cfg.spectral.lep, 'enabled',       true);
cfg.spectral.lep = defaultField(cfg.spectral.lep, 'window_sec',    [-0.1, 1.0]);
cfg.spectral.lep = defaultField(cfg.spectral.lep, 'n2_window_sec', [0.10, 0.30]);
cfg.spectral.lep = defaultField(cfg.spectral.lep, 'p2_window_sec', [0.25, 0.45]);

% ---- Hilbert phase (NEW) ----
% spectral_core first tries 09_hilbert stage (preproc_stage09_hilbert.m).
% If absent it computes phase inline using spec_compute_phase_metric.
cfg.spectral = defaultStruct(cfg.spectral, 'phase');
cfg.spectral.phase = defaultField(cfg.spectral.phase, 'enabled',        true);
cfg.spectral.phase = defaultField(cfg.spectral.phase, 'n_permutations', 5000);

% ---- Trial-spectral QC plots ----
cfg.spectral = defaultStruct(cfg.spectral, 'trial_spectral');
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'enabled',             false);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'pre_sec',             [-1.0, -0.1]);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'post_sec',            [ 0.1,  0.8]);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'fmin_hz',             1);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'fmax_hz',             40);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'max_trials',          20);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'legend_max_channels', 16);

% ---- QC / plot ----
cfg.spectral = defaultStruct(cfg.spectral, 'qc');
cfg.spectral.qc = defaultField(cfg.spectral.qc, 'plot_mode',           "summary");
cfg.spectral.qc = defaultField(cfg.spectral.qc, 'save_heatmaps',       true);
cfg.spectral.qc = defaultField(cfg.spectral.qc, 'legend_max_channels', 20);
cfg.spectral.qc = defaultField(cfg.spectral.qc, 'max_debug_trials',    5);

% ---- FOOOF bridge ----
cfg.spectral = defaultStruct(cfg.spectral, 'fooof');
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'enabled',           false);
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'python_exe',        "python3");
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'script_path',       "");
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'fmin_hz',           2);
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'fmax_hz',           40);
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'aperiodic_mode',    "fixed");
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'peak_width_limits', [1 12]);
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'max_n_peaks',       6);
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'min_peak_height',   0.0);
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'alpha_band_hz',    [8 12]);
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'verbose',           false);

if cfg.spectral.fooof.enabled
    if strlength(string(cfg.spectral.fooof.script_path)) == 0
        error('spectral_default:FooofMissingScript', ...
            'cfg.spectral.fooof.enabled=true but cfg.spectral.fooof.script_path is empty.');
    end
end

% ------------------
% Print run header
% ------------------
fprintf('[%s] SPECTRAL subjects (%d): %s\n', ...
    string(cfg.exp.id), numel(cfg.exp.subjects), mat2str(cfg.exp.subjects(:)'));
fprintf('  input_stage : %s\n', string(cfg.spectral.input_stage));
fprintf('  windows     : pre=[%.2f %.2f]s  post=[%.2f %.2f]s\n', ...
    cfg.spectral.windows.pre_sec(1),  cfg.spectral.windows.pre_sec(2), ...
    cfg.spectral.windows.post_sec(1), cfg.spectral.windows.post_sec(2));
fprintf('  LEP         : enabled=%d\n', logical(cfg.spectral.lep.enabled));
fprintf('  Phase       : enabled=%d  nPerm=%d\n', ...
    logical(cfg.spectral.phase.enabled), cfg.spectral.phase.n_permutations);

% ----------
% Dispatch
% ----------
if cfg.spectral.enabled
    spectral_core(P, cfg);
else
    fprintf('[%s] cfg.spectral.enabled=false (skipping spectral_core)\n', string(cfg.exp.id));
end
end