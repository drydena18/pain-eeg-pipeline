function preproc_stage09_hilbert(P, cfg)
% PREPROC_STAGE09_HILBERT  Standalone re-run utility for the slow-alpha Hilbert stage
% V 1.1.0
%
% Under normal operation this stage runs automatically at the end of the
% preproc_core subject loop (V2.1+).  This standalone exists so you can
% regenerate stage 09 for a subset of subjects without rerunning the full
% preprocessing pipeline, e.g.:
%
%   cfg = load_cfg(cfg_path);
%   P   = config_paths(exp_id, cfg);
%   cfg.exp.subjects = [3 7 12];   % override subject list before calling
%   preproc_stage09_hilbert(P, cfg);
%
% The output is identical to what preproc_core produces: a single-precision
% .mat file under <subRoot>/09_hilbert/sub-XXX_hilbert_phase.mat containing:
%   phase_slow    : single [nChan x nTime x nTrials]   (radians, -pi to pi)
%   t_ms          : double [1 x nTime]
%   stimOnsetIdx  : integer  (sample index of t = 0 ms)
%   chan_labels   : {nChan x 1}
%   fs            : sampling rate (Hz)
%   slow_hz       : [1 x 2]  bandpass used
%   subjid_save   : integer
%
% slow_hz resolution (mirrors preproc_core):
%   1. cfg.spectral.alpha.slow_hz  -- preferred; keeps spectral config as the
%      single source of truth for alpha band definitions
%   2. cfg.preproc.hilbert.slow_hz -- fallback; useful if spectral block is
%      absent from the JSON (default [8 10] set by normalize_preproc_defaults)
%   3. Hard fallback [8 10]        -- always safe
%
% Resume: if the output .mat already exists for a subject the subject is
% skipped, consistent with all other preprocessing stages.
%
% Requires:
%   - EEGLAB on path  (pop_eegfiltnew, eeg_checkset, pop_loadset)
%   - preproc pipeline helpers: ensure_dir, logmsg, safeClose,
%     spec_find_latest_set  (from spectral helpers, already on path when
%     called after exp01_preproc / exp01_spectral setup)

mustHave(cfg, 'exp', 'Missing cfg.exp in JSON.');
% cfg.spectral is NOT required; slow_hz is resolved safely below.

subs = cfg.exp.subjects(:);

% ---------------------------------------------------------------
% Resolve slow_hz — three-level fallback, mirrors preproc_core
% ---------------------------------------------------------------
slowHz = [8 10];   % hard fallback

if isfield(cfg, 'spectral') && isfield(cfg.spectral, 'alpha') && ...
   isfield(cfg.spectral.alpha, 'slow_hz') && ~isempty(cfg.spectral.alpha.slow_hz)
    slowHz = cfg.spectral.alpha.slow_hz;
elseif isfield(cfg, 'preproc') && isfield(cfg.preproc, 'hilbert') && ...
       isfield(cfg.preproc.hilbert, 'slow_hz') && ~isempty(cfg.preproc.hilbert.slow_hz)
    slowHz = cfg.preproc.hilbert.slow_hz;
end

% ---------------------------------------------------------------
% Resolve input stage
% ---------------------------------------------------------------
inStage = "08_base";

if isfield(cfg, 'preproc') && isfield(cfg.preproc, 'hilbert') && ...
   isfield(cfg.preproc.hilbert, 'input_stage') && ...
   strlength(string(cfg.preproc.hilbert.input_stage)) > 0
    inStage = string(cfg.preproc.hilbert.input_stage);
elseif isfield(cfg, 'spectral') && isfield(cfg.spectral, 'input_stage') && ...
       strlength(string(cfg.spectral.input_stage)) > 0
    inStage = string(cfg.spectral.input_stage);
end

fprintf('[STAGE09] Standalone Hilbert | subjects=%d  band=[%.0f %.0f] Hz  input=%s\n', ...
    numel(subs), slowHz(1), slowHz(2), inStage);

for i = 1:numel(subs)
    subjid = subs(i);

    subRoot   = fullfile(string(P.RUN_ROOT), sprintf('sub-%03d', subjid));
    inDir     = fullfile(subRoot, char(inStage));
    outDir    = fullfile(subRoot, '09_hilbert');
    outLogDir = fullfile(outDir, 'logs');

    ensure_dir(outDir);
    ensure_dir(outLogDir);

    % Use preproc-style logging (fopen/logmsg/safeClose) to avoid
    % dependency on spectral helpers being on the path
    logPath = fullfile(outLogDir, sprintf('sub-%03d_stage09_hilbert.log', subjid));
    logf = fopen(logPath, 'w');
    if logf < 0
        warning('preproc_stage09_hilbert:LogFail', 'Could not open log: %s', logPath);
        logf = 1;
    end
    cobj = onCleanup(@() safeClose(logf));   

    logmsg(logf, '===== STAGE09_HILBERT (standalone) START sub-%03d =====', subjid);
    logmsg(logf, 'band=[%.0f %.0f] Hz  input=%s', slowHz(1), slowHz(2), inStage);

    % ---------------------------------------------------------------
    % Resume guard — skip if output already exists
    % ---------------------------------------------------------------
    matPath = fullfile(outDir, sprintf('sub-%03d_hilbert_phase.mat', subjid));

    if exist(matPath, 'file')
        logmsg(logf, '[SKIP] Output already exists (resume): %s', matPath);
        logmsg(logf, '===== STAGE09_HILBERT DONE sub-%03d =====', subjid);
        continue;
    end

    % ---------------------------------------------------------------
    % Load 08_base
    % ---------------------------------------------------------------
    if ~exist(inDir, 'dir')
        logmsg(logf, '[WARN] Input dir missing: %s (skipping)', inDir);
        continue;
    end

    inSet = spec_find_latest_set(inDir, cfg.exp.out_prefix, subjid);
    if strlength(inSet) == 0
        logmsg(logf, '[WARN] No .set found in %s (skipping)', inDir);
        continue;
    end

    logmsg(logf, '[LOAD] %s', inSet);
    inSet = char(inSet);
    [inFolder, inName, inExt] = fileparts(inSet);
    EEG = pop_loadset('filename', [inName inExt], 'filepath', inFolder);
    EEG = eeg_checkset(EEG);

    if EEG.trials <= 1
        logmsg(logf, '[WARN] EEG not epoched (trials=%d). Stage 09 requires epochs. Skipping.', EEG.trials);
        continue;
    end

    fs    = EEG.srate;
    nChan = EEG.nbchan;
    nTr   = EEG.trials;
    nTime = EEG.pnts;

    % Stimulus-onset sample (t closest to 0 ms)
    timesSec   = double(EEG.times(:))' / 1000;
    [~, t0idx] = min(abs(timesSec));
    logmsg(logf, '[STAGE09] stimOnsetIdx=%d  t=%.3f ms', t0idx, EEG.times(t0idx));

    % ---------------------------------------------------------------
    % Zero-phase FIR bandpass
    % ---------------------------------------------------------------
    bw       = slowHz(2) - slowHz(1);
    minOrder = 3 * round(fs / bw);
    logmsg(logf, '[STAGE09] Bandpass [%.0f %.0f] Hz  FIR order >= %d', ...
        slowHz(1), slowHz(2), minOrder);

    try
        EEGfilt = pop_eegfiltnew(EEG, slowHz(1), slowHz(2));
        EEGfilt = eeg_checkset(EEGfilt);
    catch ME
        logmsg(logf, '[WARN] FIR bandpass failed: %s (skipping sub-%03d)', ME.message, subjid);
        continue;
    end

    % ---------------------------------------------------------------
    % Hilbert transform per channel x trial
    % ---------------------------------------------------------------
    phaseSlow = zeros(nChan, nTime, nTr, 'single');

    logmsg(logf, '[STAGE09] Computing Hilbert phase (%d chans x %d trials)...', nChan, nTr);

    for t = 1:nTr
        X = double(EEGfilt.data(:, :, t));   % [nChan x nTime]
        for ch = 1:nChan
            z = hilbert(X(ch, :)');           % analytic signal [nTime x 1]
            phaseSlow(ch, :, t) = single(angle(z));
        end
    end

    logmsg(logf, '[STAGE09] Hilbert done.  phase_slow size: %s', mat2str(size(phaseSlow)));

    % ---------------------------------------------------------------
    % Save
    % ---------------------------------------------------------------
    phase_slow   = phaseSlow;              %#ok<NASGU>
    t_ms         = EEG.times;             %#ok<NASGU>
    stimOnsetIdx = t0idx;                 %#ok<NASGU>
    chan_labels  = {EEG.chanlocs.labels}; %#ok<NASGU>
    slow_hz      = slowHz;                %#ok<NASGU>
    subjid_save  = subjid;                %#ok<NASGU>

    save(matPath, 'phase_slow', 't_ms', 'stimOnsetIdx', ...
        'chan_labels', 'slow_hz', 'subjid_save', 'fs', '-v7.3');

    logmsg(logf, '[STAGE09] Saved -> %s', matPath);
    logmsg(logf, '===== STAGE09_HILBERT DONE sub-%03d =====', subjid);
end

fprintf('[STAGE09] Standalone complete.\n');
end