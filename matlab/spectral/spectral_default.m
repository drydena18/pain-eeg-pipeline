function cfg = spectral_default(P, cfg_in, subjects_override)
% SPECTRAL_DEFAULT Validate/normalise cfg then dispatch to spectral_core
% V 2.1.0
%
% V2.1.0: features.windows.whole_sec (default [-1 2] s) replaces 'whole =
% full epoch' now that epochs run -2.2 to 2.0 s; warns when a window lies
% within half a band-power filter length of an epoch edge.
%
% V2.0.0 changes vs V1.1.0 (channel/source parity):
%   - Feature parameters come from the top-level cfg.features block, the
%     single source of truth shared with python/source/source_default.py:
%     windows, bands, band-power filter, Welch PSD, ERD guard, LEP windows,
%     phase reference and FOOOF. They are normalised into cfg.spectral.feat,
%     which is all spectral_core reads for features. Legacy cfg.spectral.*
%     values that conflict with cfg.features are reported and ignored.
%   - Features are validated against the PREPROCESSING baseline in the same
%     JSON (identical checks to source_default.py): windows inside the epoch,
%     bands inside the pass band, band-power FIR shorter than the epoch.
%   - cfg.spectral.fooof keeps only python_exe / script_path / enabled.
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

% ---------------------------------------------
% Shared features (single source of truth)
% ---------------------------------------------
userFeat = struct();
if isfield(cfg, 'features') && isstruct(cfg.features)
    userFeat = cfg.features;
end
F = local_merge(local_feature_defaults(), userFeat);
F = local_shape(F);
local_validate_features(F, cfg);
local_report_conflicts(cfg.spectral, F);
cfg.features = F;

% ---------------------------
% Spectral block (non-feature settings)
% ---------------------------
cfg.spectral = defaultField(cfg.spectral, 'enabled',     true);
cfg.spectral = defaultField(cfg.spectral, 'input_stage', "08_base");
cfg.spectral.feat = F;

% LEP waveform storage span (peaks use F.lep)
cfg.spectral = defaultStruct(cfg.spectral, 'lep');
cfg.spectral.lep = defaultField(cfg.spectral.lep, 'enabled',    true);
cfg.spectral.lep = defaultField(cfg.spectral.lep, 'window_sec', [-0.1, 1.0]);

% Phase (MATLAB-only extras: r_cl permutation test against ratings)
cfg.spectral = defaultStruct(cfg.spectral, 'phase');
cfg.spectral.phase = defaultField(cfg.spectral.phase, 'enabled',        true);
cfg.spectral.phase = defaultField(cfg.spectral.phase, 'n_permutations', 5000);

% QC plots
cfg.spectral = defaultStruct(cfg.spectral, 'trial_spectral');
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'enabled',             false);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'pre_sec',             F.windows.pre_sec);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'post_sec',            F.windows.post_sec);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'fmin_hz',             F.psd.fmin_hz);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'fmax_hz',             F.psd.fmax_hz);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'max_trials',          20);
cfg.spectral.trial_spectral = defaultField(cfg.spectral.trial_spectral, 'legend_max_channels', 16);

cfg.spectral = defaultStruct(cfg.spectral, 'qc');
cfg.spectral.qc = defaultField(cfg.spectral.qc, 'plot_mode',           "summary");
cfg.spectral.qc = defaultField(cfg.spectral.qc, 'save_heatmaps',       true);
cfg.spectral.qc = defaultField(cfg.spectral.qc, 'legend_max_channels', 20);
cfg.spectral.qc = defaultField(cfg.spectral.qc, 'max_debug_trials',    5);

% Legacy keys some QC plotting helpers still read
cfg.spectral.alpha = struct('alpha_hz', F.bands.alpha_hz, 'slow_hz', F.bands.slow_hz, 'fast_hz', F.bands.fast_hz);
cfg.spectral.windows = F.windows;

% FOOOF environment (parameters are in F.fooof)
cfg.spectral = defaultStruct(cfg.spectral, 'fooof');
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'python_exe',  "python3");
cfg.spectral.fooof = defaultField(cfg.spectral.fooof, 'script_path', "");
if F.fooof.enabled && strlength(string(cfg.spectral.fooof.script_path)) == 0
    error('spectral_default:FooofMissingScript', ...
        'features.fooof.enabled = true but cfg.spectral.fooof.script_path is empty.');
end

% ------------------
% Print run header
% ------------------
fprintf('[%s] SPECTRAL subjects (%d): %s\n', ...
    string(cfg.exp.id), numel(cfg.exp.subjects), mat2str(cfg.exp.subjects(:)'));
fprintf('  input_stage : %s\n', string(cfg.spectral.input_stage));
fprintf('  windows     : whole=[%.2f %.2f]s  pre=[%.2f %.2f]s  post=[%.2f %.2f]s\n', ...
    F.windows.whole_sec(1), F.windows.whole_sec(2), F.windows.pre_sec(1), F.windows.pre_sec(2), F.windows.post_sec(1), F.windows.post_sec(2));
fprintf('  band power  : filter-Hilbert, slow %s / fast %s Hz, trans_bw %.2f Hz\n', ...
    mat2str(F.bands.slow_hz), mat2str(F.bands.fast_hz), F.band_power.trans_bw_hz);
fprintf('  LEP         : N2 %s s, P2 %s s\n', mat2str(F.lep.n2_window_sec), mat2str(F.lep.p2_window_sec));
fprintf('  FOOOF       : enabled=%d (subject x channel trial-mean PSD)\n', logical(F.fooof.enabled));

% ----------
% Dispatch
% ----------
if cfg.spectral.enabled
    spectral_core(P, cfg);
else
    fprintf('[%s] cfg.spectral.enabled=false (skipping spectral_core)\n', string(cfg.exp.id));
end
end

% ==========================================================================
% Local: shared feature defaults (mirror of source_default.FEATURE_DEFAULTS)
% ==========================================================================
function F = local_feature_defaults()
F = struct();
F.windows    = struct('whole_sec', [-1.0, 2.0], 'pre_sec', [-1.0, -0.1], 'post_sec', [0.1, 0.8]);
F.bands      = struct('alpha_hz', [8, 12], 'slow_hz', [8, 10], 'fast_hz', [10, 12]);
F.band_power = struct('method', 'filter_hilbert', 'trans_bw_hz', 1.5);
F.psd        = struct('fmin_hz', 1, 'fmax_hz', 40, 'window_sec', 2.0, 'overlap_frac', 0.5, 'df_target_hz', 0.25);
F.erd        = struct('eps0_frac', 1e-3, 'p5_percentile', 5);
F.lep        = struct('n2_window_sec', [0.15, 0.35], 'p2_window_sec', [0.25, 0.50]);
F.phase      = struct('post_ref_sec', 0.2);
F.fooof      = struct('enabled', true, 'aperiodic_mode', 'fixed', 'peak_width_limits', [1, 12], ...
                      'max_n_peaks', 6, 'min_peak_height', 0.1, 'peak_threshold', 2.0, 'freq_range', [1, 40]);
end

function out = local_merge(out, user)
fn = fieldnames(user);
for i = 1:numel(fn)
    k = fn{i};
    if startsWith(k, 'x_') || startsWith(k, '_')   % "_comment" keys -> x_comment after jsondecode
        continue;
    end
    if isstruct(user.(k)) && isfield(out, k) && isstruct(out.(k))
        out.(k) = local_merge(out.(k), user.(k));
    else
        out.(k) = user.(k);
    end
end
end

function F = local_shape(F)
% jsondecode returns column vectors; make every 2-element range a row.
groups = fieldnames(F);
for g = 1:numel(groups)
    if ~isstruct(F.(groups{g})), continue; end
    fn = fieldnames(F.(groups{g}));
    for i = 1:numel(fn)
        v = F.(groups{g}).(fn{i});
        if isnumeric(v) && numel(v) > 1
            F.(groups{g}).(fn{i}) = double(v(:)');
        end
    end
end
F.fooof.enabled = logical(F.fooof.enabled);
F.fooof.aperiodic_mode = char(string(F.fooof.aperiodic_mode));
end

function local_validate_features(F, cfg)
% Mirror of source_default.validate_features_against_preproc.
tmin = -1.0; tmax = 2.0; hp = 0.5; lp = 40; fs = 500;
if isfield(cfg, 'preproc')
    pp = cfg.preproc;
    if isfield(pp, 'epoch') && isfield(pp.epoch, 'tmin_sec'), tmin = pp.epoch.tmin_sec; end
    if isfield(pp, 'epoch') && isfield(pp.epoch, 'tmax_sec'), tmax = pp.epoch.tmax_sec; end
    if isfield(pp, 'filter') && isfield(pp.filter, 'highpass_hz'), hp = pp.filter.highpass_hz; end
    if isfield(pp, 'filter') && isfield(pp.filter, 'lowpass_hz'),  lp = pp.filter.lowpass_hz;  end
    if isfield(pp, 'resample') && isfield(pp.resample, 'target_hz') && ~isempty(pp.resample.target_hz)
        fs = pp.resample.target_hz;
    end
end

w = {'whole_sec', 'pre_sec', 'post_sec'};
for i = 1:3
    r = F.windows.(w{i});
    if ~(tmin <= r(1) && r(1) < r(2) && r(2) <= tmax)
        error('spectral_default:WindowOutsideEpoch', ...
            'features.windows.%s = %s is not inside the epoch [%g %g] s.', w{i}, mat2str(r), tmin, tmax);
    end
end

b = F.bands;
if ~(b.slow_hz(2) == b.fast_hz(1) && b.slow_hz(1) == b.alpha_hz(1) && b.fast_hz(2) == b.alpha_hz(2))
    warning('spectral_default:BandsDoNotTile', 'features.bands: slow + fast do not tile alpha.');
end
bn = {'alpha_hz', 'slow_hz', 'fast_hz'};
for i = 1:3
    r = b.(bn{i});
    if ~(hp < r(1) && r(1) < r(2) && r(2) < lp)
        error('spectral_default:BandOutsidePassband', ...
            'features.bands.%s = %s is outside the preprocessing pass band (%g-%g Hz).', bn{i}, mat2str(r), hp, lp);
    end
end

nTaps = floor(3.3 * fs / F.band_power.trans_bw_hz + 0.5);
if mod(nTaps, 2) == 0, nTaps = nTaps + 1; end
nEpoch = round((tmax - tmin) * fs) + 1;
if nTaps > nEpoch
    error('spectral_default:FilterTooLong', ...
        'Band-power FIR (%d taps at %g Hz) is longer than the epoch (%d samples).', nTaps, fs, nEpoch);
end
% Edge zone: samples within half a filter length of an epoch edge depend on
% the mirror padding (pre-window power ~9 % low when the window touches the edge).
half = (nTaps - 1) / 2 / fs;
for i = 1:3
    r = F.windows.(w{i});
    if r(1) - tmin < half - 1e-9 || tmax - r(2) < half - 1e-9
        warning('spectral_default:WindowInEdgeZone', ...
            'features.windows.%s = %s lies within %.2f s of the epoch edge [%g %g] s; band power there depends on filter padding.', ...
            w{i}, mat2str(r), half, tmin, tmax);
    end
end

if F.psd.fmax_hz > lp
    warning('spectral_default:PsdAboveLowpass', 'features.psd.fmax_hz = %g exceeds the %g Hz low-pass.', F.psd.fmax_hz, lp);
end

post = F.windows.post_sec;
ln = {'n2_window_sec', 'p2_window_sec'};
for i = 1:2
    r = F.lep.(ln{i});
    if ~(post(1) <= r(1) && r(1) < r(2) && r(2) <= post(2))
        warning('spectral_default:LepOutsidePost', 'features.lep.%s = %s is not inside the post window %s.', ...
            ln{i}, mat2str(r), mat2str(post));
    end
end

if ~(tmin < F.phase.post_ref_sec && F.phase.post_ref_sec < tmax)
    error('spectral_default:PhaseRefOutsideEpoch', 'features.phase.post_ref_sec is outside the epoch.');
end
end

function local_report_conflicts(S, F)
% Warn when legacy cfg.spectral.* values disagree with cfg.features.
pairs = {
    'windows', 'pre_sec',  F.windows.pre_sec
    'windows', 'post_sec', F.windows.post_sec
    'alpha',   'alpha_hz', F.bands.alpha_hz
    'alpha',   'slow_hz',  F.bands.slow_hz
    'alpha',   'fast_hz',  F.bands.fast_hz
    'lep',     'n2_window_sec', F.lep.n2_window_sec
    'lep',     'p2_window_sec', F.lep.p2_window_sec
};
for i = 1:size(pairs, 1)
    g = pairs{i, 1}; k = pairs{i, 2};
    if isfield(S, g) && isstruct(S.(g)) && isfield(S.(g), k)
        v = double(S.(g).(k)(:)');
        if ~isequal(v, pairs{i, 3})
            warning('spectral_default:LegacyConflict', ...
                'cfg.spectral.%s.%s = %s conflicts with cfg.features (%s); using features.', ...
                g, k, mat2str(v), mat2str(pairs{i, 3}));
        end
    end
end
end
