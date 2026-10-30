function cfg = pp_defaults()
% PP_DEFAULTS Generic config for the two-mode (unattended + manual) runner.
%   Normally built by pp_from_experiment('expXX'), which fills in
%   the experiment and wires the stage hooks from pp_hooks.m.
%
%   Stage hook signature: [EEG, qc] = fn(EEG, cfg, subj)

cfg = struct();

% ---- Subjects & paths --------------------------------------------------
cfg.subjects     = {};        % e.g. {'sub-001','sub-002'}
cfg.paths.out    = '';        % experiment root; per-subject folders live under it
cfg.paths.code   = '';        % repo root; git commit hash is logged if set
cfg.paths.prefix = '';        % file prefix, e.g. '142BB_62_' -> 142BB_62_006_...
% Per-subject folders (<out>/sub-XXX/<folder>); logs at experiment level.
cfg.paths.sub = struct('ica', '06_ica_ready', 'review', '06_ica_review', ...
                       'final', '08_base', 'logs', '_pp_logs');

% Optional: [rt, cleanup] = fn(cfg, subj, mode). Each mode puts rt in
% cfg.rt before calling stage hooks (log handle, stage folders, ...).
cfg.runtime_fn = [];

% ---- Mode 1: unattended preprocessing -> ICA-ready ------------------------
cfg.mode1.version       = 'v1';
cfg.mode1.load_fn       = @(cfg, subj) error('pp:todo', 'Set cfg.mode1.load_fn.');
cfg.mode1.stages        = struct('name', {}, 'fn', {});
cfg.mode1.clean_copy_fn = @(EEG, cfg, subj) error('pp:todo', 'Set cfg.mode1.clean_copy_fn.');

% Rank bookkeeping - must describe the pipeline as it actually runs.
cfg.mode1.avg_ref             = true;
cfg.mode1.reref_before_interp = true;
cfg.mode1.rank_eig_tol        = 1e-7;

cfg.mode1.ica_extended = 1;            % extended Infomax
cfg.mode1.ica_seed     = 'subject';    % rng(<subject number>) like preproc_core, or a number

% ---- Mode 2a: automated review prep ---------------------------------------
cfg.mode2a.version         = 'v1';
cfg.mode2a.prep_fn         = [];       % [EEG, suggested, qc] = fn(EEG, cfg, subj); [] = built-in
cfg.mode2a.iclabel_version = 'default';
cfg.mode2a.n_show          = 50;
cfg.mode2a.flag = struct('Muscle', 0.9, 'Eye', 0.9, 'Heart', 0.9, ...
                         'LineNoise', 0.9, 'ChannelNoise', 0.9);   % built-in prep only

% ---- Manual review --------------------------------------------------------
cfg.review.reviewer   = char(java.lang.System.getProperty('user.name'));
cfg.review.ui_fn      = @pp_review_ui_default;
cfg.review.autorun_2b = true;

% ---- Mode 2b: apply decisions -> final dataset ----------------------------
cfg.mode2b.version = 'v1';
cfg.mode2b.stages  = struct('name', {}, 'fn', {});
end