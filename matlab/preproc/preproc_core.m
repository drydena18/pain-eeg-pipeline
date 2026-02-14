function preproc_core(P, cfg)
% PREPROC_CORE Execute preprocessing according to cfg (from JSON)
<<<<<<< HEAD
% V 2.0.0
%
% Folder scheme (per subject):
%   PROJ_ROOT/<exp_out>/sub-XXX/01_filter
%   PROJ_ROOT/<exp_out>/sub-XXX/02_notch
%   PROJ_ROOT/<exp_out>/sub-XXX/03_resample
%   PROJ_ROOT/<exp_out>/sub-XXX/04_reref
%   PROJ_ROOT/<exp_out>/sub-XXX/05_initrej
%   PROJ_ROOT/<exp_out>/sub-XXX/06_ica
%   PROJ_ROOT/<exp_out>/sub-XXX/07_epoch
%   PROJ_ROOT/<exp_out>/sub-XXX/08_base
%   PROJ_ROOT/<exp_out>/sub-XXX/LOGS    (csv/txt/png + .log)
%   PROJ_ROOT/<exp_out>/sub-XXX/QC      (QC packets / structured QC
%   outputs)
%
% Requires:
%   - EEGLAB on path
%   - BIOSIG plugin for pop_biosig (for .eeg/.bdf)
%   - ICLabel plugin if enabled

subs = cfg.exp.subjects(:);

% Resolve experiment output folder name (stable)
expRoot = string(P.RUN_ROOT);
ensure_dir(expRoot);

for i = 1:numel(subs)
    subjid = subs(i);
    tags = {}; % cumulative tags for this subject

    % --------------------------
    % Per-subject folder scheme
    % --------------------------
    subRoot = fullfile(expRoot, sprintf('sub-%03d', subjid));

    ST = struct();
    ST.FILTER   = fullfile(subRoot, '01_filter');
    ST.NOTCH    = fullfile(subRoot, '02_notch');
    ST.RESAMPLE = fullfile(subRoot, '03_resample');
    ST.REREF    = fullfile(subRoot, '04_reref');
    ST.INITREJ  = fullfile(subRoot, '05_initrej');
    ST.ICA      = fullfile(subRoot, '06_ica');
    ST.EPOCH    = fullfile(subRoot, '07_epoch');
    ST.BASE     = fullfile(subRoot, '08_base');
    
    LOGS = fullfile(subRoot, 'LOGS');
    QC   = fullfile(subRoot, 'QC');

    ensure_dir(subRoot);
    ensure_dir(LOGS);
    ensure_dir(QC);

    fns = fieldnames(ST);
    for k = 1:numel(fns)
        ensure_dir(ST.(fns{k}));
    end

    % ----------------------
    % Logging (per subject)
    % ----------------------
    logPath = fullfile(LOGS, sprintf('sub-%03d_preproc.log', subjid));
    logf = fopen(logPath, 'w');
    if logf < 0
        warning('preproc_core:LogOpenFail', 'Could not open log file: %s', logPath);
        logf = 1; % fallback to stdout
    end
    cleanupObj = onCleanup(@() safeClose(logf));

    logmsg(logf, '===== PREPROC START sub-%03d =====', subjid);
    logmsg(logf, 'Experiment: %s', string(cfg.exp.id));
    logmsg(logf, 'Out Prefix: %s', string(cfg.exp.out_prefix));
    logmsg(logf, 'Output Root: %s', subRoot);

    % -----------------------
    % Resolve raw input file
    % -----------------------
    rawPath = resolve_raw_file(P, cfg, subjid);
    if isempty(rawPath) || ~exist(rawPath, 'file')
        logmsg(logf, '[WARN] Raw file not found for sub-%03d. Skipping.', subjid);
        continue;
    end
    logmsg(logf, 'Raw File: %s', rawPath);

    % -------------------------------------------------
    % Import raw (.eeg/.bdf) -> EEG struct (no saving)
    % -------------------------------------------------
=======
% V 1.2.0
%
% Pipeline:
%   import -> filter -> notch -> resample -> reref -> initrej (manual
%   interp) -> ica (trained on optional clean copy) -> iclabel (suggest) ->
%   manual ic reject -> epoch -> baseline
%
% Naming:
%   base: 26BB_64_001.set (conceptual; not saved)
%   storage saves: use cumulative tags, e.g.,
%   26BB_64_001_fir_notch60_rs500_reref.set
%
% Requirements:
%   - EEGLAB on path
%   - BIOSIG plugin for pop_biosig (for .eeg)
%   - ICLabel plugin if cfg.preproc.ica.iclabel.enabled = true

subs = cfg.exp.subjects(:);

for i = 1:numel(subs)
    subjid = subs(i);
    tags = {}; % cumulative tags for this subject

    % ---------------------
    % Logging (per subject)
    % ---------------------
    logPath = fullfile(P.STAGE.LOGS, sprintf('sub-%03d_preproc.log', subjid));
    logf = fopen(logPath, 'w');
    if logf < 0
        warning('preproc_core:LogOpenFail', 'Could not open log file: %s', logPath);
        logf = 1; % fallback to stdout
    end
    cleanupObj = onCleanup(@() safeClose(logf)); %#ok<NASGU>

    logmsg(logf, '===== PREPROC START sub-%03d =====', subjid);
    logmsg(logf, 'Experiment: %s', string(cfg.exp.id));
    logmsg(logf, 'Out prefix: %s', string(cfg.exp.out_prefix));

    % ----------------------
    % Resolve raw input file
    % ----------------------
    rawPath = resolve_raw_file(P, cfg, subjid);
    if isempty(rawPath) || ~exist(rawPath, 'file')
        logmsg(logf, '[WARN] Raw file not found for sub-%03d. Skipping.', subjid);
        continue;
    end
    logmsg(logf, 'Raw file: %s', rawPath);

    % ------------------------------------------------
    % Inport raw (.eeg) -> EEG struct (no saving here)
    % ------------------------------------------------
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
    EEG = pop_biosig(rawPath);
    EEG = eeg_checkset(EEG);

    EEG = normalize_chan_labels(EEG);

<<<<<<< HEAD
    % -------------------
    % Montage (optional)
    % -------------------
    if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'enabled') && cfg.exp.montage.enabled
        EEG = apply_montage_biosemi_from_csv(P, cfg, EEG, logf, subjid, LOGS);
        EEG = eeg_checkset(EEG);
    end

    % --------------------------------------
    % Channel Locations (.elp) if requested
    % --------------------------------------
    if isfield(cfg.exp, 'channel_locs') && isfield(cfg.exp.channel_locs, 'use_elp') && cfg.exp.channel_locs.use_elp
        elpPath = "";
        if isfield(P, 'CORE') && isfield(P.CORE, 'ELP_FILE')
            elpPath = string(P.CORE.ELP_FILE);
        end

        if strlength(elpPath) > 0 && exist(elpPath, 'file')
            try
                EEG = pop_chanedit(EEG, 'lookup', char(elpPath));
                EEG = eeg_checkset(EEG);
                logmsg(logf, 'Loaded chanlocs from ELP: %s', elpPath);
            catch ME
                logmsg(logf, '[WARN] Failed to load chanlocs: %s', ME.message);
            end
        else
            logmsg(logf, '[WARN] ELP file missing: %s', elpPath);
        end
=======
    % Apply experiment-specific montage (e.g., Biosemi A/B -> 10-20)
    if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'enabled') && cfg.exp.montage.enabled
        EEG = apply_montage_biosemi_from_csv(P, cfg, EEG, logf, subjid);
        EEG = eeg_checkset(EEG);
    end

    % -------------------------------------
    % Channel locations (.elp) if requested
    % -------------------------------------
    if isfield(cfg.exp, 'channel_locs') && isfield(cfg.exp.channel_locs, 'use_elp') && cfg.exp.channel_locs.use_elp
        if isfield(P, 'CORE') && isfield(P.CORE, 'ELP_FILE') && exist(P.CORE.ELP_FILE, 'file')
            try
                EEG = pop_chanedit(EEG, 'lookup', P.CORE.ELP_FILE);
                EEG = eeg_checkset(EEG);
                logmsg(logf, 'Loaded chanlocs from ELP: %s', P.CORE.ELP_FILE);
            catch ME
                logmsg(logf, '[WARN] Failed to load chanlocs: %s', ME.message);
            end
        else
            logmsg(logf, '[WARN] ELP file missing: %s', P.CORE.ELP_FILE);
        end
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
    end

    % Deterministic seed per subject (repro ICA runs)
    rng(double(subjid), 'twister');
    logmsg(logf, '[RNG] rng(%d, twister)', subjid);

<<<<<<< HEAD
    %% ------
    % FILTER
    % -------
    if cfg.preproc.filter.enabled
        nextTag = char(string(cfg.preproc.filter.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.FILTER, P, subjid, tags, nextTag, logf);
=======
    % --------------
    % FILTER (HP/LP)
    % --------------
    if cfg.preproc.filter.enabled
        nextTag = char(string(cfg.preproc.filter.tag));
        [EEG, tags, didLoad] = maybe_load_stage(P.STAGE.FILTER, P, subjid, tags, nextTag, logf);
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
        if ~didLoad
            hp = cfg.preproc.filter.highpass_hz;
            lp = cfg.preproc.filter.lowpass_hz;
            logmsg(logf, '[FILTER] %s hp = %.3f lp = %.3f', string(cfg.preproc.filter.type), hp, lp);

            EEG = pop_eegfiltnew(EEG, hp, lp);
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
<<<<<<< HEAD
            save_stage(ST.FILTER, P, subjid, tags, EEG, logf);
        end
    end

    %% ----------------------------
    % NOTCH (bandstop via revfilt)
    % -----------------------------
    if cfg.preproc.notch.enabled
        nextTag = char(string(cfg.preproc.notch.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.NOTCH, P, subjid, tags, nextTag, logf);
=======
            save_stage(P.STAGE.FILTER, P, subjid, tags, EEG, logf);
        end
    end

    % ----------------------------
    % NOTCH (bandstop via revfilt)
    % ----------------------------
    if cfg.preproc.notch.enabled
        nextTag = char(string(cfg.preproc.notch.tag));
        [EEG, tags, didLoad] = maybe_load_stage(P.STAGE.FILTER, P, subjid, tags, nextTag, logf);
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
        if ~didLoad
            f0 = cfg.preproc.notch.freq_hz;
            bw = cfg.preproc.notch.bw_hz;
            logmsg(logf, '[NOTCH] f0 = %g bw = %g', f0, bw);

<<<<<<< HEAD
            EEG = pop_eegfiltnew(EEG, f0-bw, f0+bw, [], 1); % revfilt = 1 -> bandstop
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
            save_stage(ST.NOTCH, P, subjid, tags, EEG, logf);
        end
    end

    %% --------
    % RESAMPLE
    % ---------
    if cfg.preproc.resample.enabled
        nextTag = char(string(cfg.preproc.resample.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.RESAMPLE, P, subjid, tags, nextTag, logf);
=======
            EEG = pop_eegfiltnew(EEG, f0-bw, f0+bw, [], 1); % revfilt = 1 => bandstop
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
            save_stage(P.STAGE.FILTER, P, subjid, tags, EEG, logf);
        end
    end

    % --------
    % RESAMPLE
    % --------
    if cfg.preproc.resample.enabled && ~isempty(cfg.preproc.resample.target_hz)
        nextTag = char(string(cfg.preproc.resample.tag));
        [EEG, tags, didLoad] = maybe_load_stage(P.STAGE.RESAMPLE, P, subjid, tags, nextTag, logf);
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
        if ~didLoad
            targetFs = cfg.preproc.resample.target_hz;
            if EEG.srate ~= targetFs
                logmsg(logf, '[RESAMPLE] %g -> %g Hz', EEG.srate, targetFs);
                EEG = pop_resample(EEG, targetFs);
                EEG = eeg_checkset(EEG);

                tags{end+1} = nextTag;
<<<<<<< HEAD
                save_stage(ST.RESAMPLE, P, subjid, tags, EEG, logf);
=======
                save_stage(P.STAGE.RESAMPLE, P, subjid, tags, EEG, logf);
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
            else
                logmsg(logf, '[RESAMPLE] Already at %g Hz (skipping).', targetFs);
            end
        end
    end

<<<<<<< HEAD
    %% -----------
    % REREFERENCE
    % ------------
    if cfg.preproc.reref.enabled
        nextTag = char(string(cfg.preproc.reref.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.REREF, P, subjid, tags, nextTag, logf);
=======
    % -----------
    % REREFERENCE
    % -----------
    if cfg.preproc.reref.enabled
        nextTag = char(string(cfg.preproc.reref.tag));
        [EEG, tags, didLoad] = maybe_load_stage(P.STAGE.REREF, P, subjid, tags, nextTag, logf);
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
        if ~didLoad
            logmsg(logf, '[REREF] mode = %s', string(cfg.preproc.reref.mode));

            mode = lower(string(cfg.preproc.reref.mode));
            if mode == "average"
                EEG = pop_reref(EEG, []);
            elseif mode == "channels"
                chans = cfg.preproc.reref.channels;
                if isempty(chans)
                    error('reref.mode = "channels" but cfg.preproc.reref.channels is empty.');
                end
                EEG = pop_reref(EEG, chans);
            else
                error('Unsupported reref.mode: %s', mode);
            end
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
<<<<<<< HEAD
            save_stage(ST.REREF, P, subjid, tags, EEG, logf);
        end
    end

    %% ---------------------------------------------------
    % INITREJ (suggest bad chans + plots + manual interp)
    % ----------------------------------------------------
    if cfg.preproc.initrej.enabled
        nextTag = char(string(cfg.preproc.initrej.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.INITREJ, P, subjid, tags, nextTag, logf);
        if ~didLoad
            logmsg(logf, '[INITREJ] Suggesting bad channels + manual spherical interpolation.');

=======
            save_stage(P.STAGE.REREF, P, subjid, tags, EEG, logf);
        end
    end

    % ------------------------------------------------
    % INITREJ (manual interp + suggested bad channels + histograms + PSD metrics)
    % ------------------------------------------------
    if cfg.preproc.initrej.enabled
        nextTag = char(string(cfg.preproc.initrej.tag));
        [EEG, tags, didLoad] = maybe_load_stage(P.STAGE.INITREJ, P, subjid, tags, nextTag, logf);
        if ~didLoad
            logmsg(logf, '[INITREJ] Suggesting bad channels + manual spherical interpolation.');

            % 1) Suggest bad chans (prob + kurtosis)
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
            [badChans, reasons, metrics] = suggest_bad_channels(EEG);

            if ~isfield(EEG, 'etc') || isempty(EEG.etc)
                EEG.etc = struct();
<<<<<<< HEAD
            end

            EEG.etc.initrej = struct();
            EEG.etc.initrej.suggested_badchans = badChans;
            EEG.etc.initrej.reasons = reasons;
            EEG.etc.initrej.metrics = metrics;

            % Channel PSD metrics + CSV + topoplots
            try
                chanPSD = compute_channel_psd_metrics(EEG);
                EEG.etc.initrej.chan_psd = chanPSD;
                write_channel_psd_csv(LOGS, subjid, EEG, chanPSD);
                save_chan_psd_topos(LOGS, subjid, EEG, chanPSD);
                logmsg(logf, '[INITREJ] Channel PSD metrics saved (CSV + topos).');
            catch ME
                logmsg(logf, '[WARN] Channel PSD metrics failed: %s', ME.message);
            end

            % QC plots -> LOGS
            try
                make_initrej_plots(LOGS, subjid, EEG, metrics, badChans);
                logmsg(logf, '[INITREJ] Saved QC plots to LOGS.');
            catch ME
                logmsg(logf, '[WARN] INITREJ plotting failed: %s', ME.message);
            end

            % Print suggestions
            if isempty(badChans)
                logmsg(logf, '[INITREJ] No channel suggested.');
            else
                logmsg(logf, '[INITREJ] Suggested bad channels: %s', vec2str(badChans));
                for k = 1:numel(badChans)
                    ch = badChans(k);
                    lbl = safe_chan_label(EEG, ch);
                    logmsg(logf, '  - Ch %d (%s): %s', ch, lbl, reasons{k});
                end
            end

            % Manual decision: which to interp
            interpChans = prompt_channel_interp(EEG, badChans);

            if isempty(interpChans)
                logmsg(logf, '[INITREJ] No channels interpolated (manual decision).');
=======
            end
            EEG.etc.initrej = struct();
            EEG.etc.initrej.suggested_badchans = badChans;
            EEG.etc.initrej.reasons = reasons;
            EEG.etc.initrej.metrics = metrics;

            % Channel PSD metrics + CSV + topoplots
            try
                chanPSD = compute_channel_psd_metrics(EEG);
                EEG.etc.initrej.chan_psd = chanPSD;
                write_channel_psd_csv(P, subjid, EEG, chanPSD);
                save_chan_psd_topos(P, subjid, EEG, chanPSD);
                logmsg(logf, '[INITREJ] Channel PSD metrics saved (CSV + topo is chanlocs).');
            catch ME
                logmsg(logf, '[WARN] Channel PSD metrics failed: %s', ME.message);
            end

            % 2) QC plots (hist/bar/topo/PSD) saved to logs
            try
                make_initrej_plots(P, subjid, EEG, metrics, badChans);
                logmsg(logf, '[INITREJ] Saved QC plots to logs.\n');
            catch ME
                logmsg(logf, '[WARN] INITREJ plotting failed: %s', ME.message);
            end

            % 3) Print suggestions
            if isempty(badChans)
                logmsg(logf, '[INITREJ] No channels suggested.');
            else
                logmsg(logf, '[INITREJ] Suggested bad channels: %s', vec2str(badChans));
                for k = 1:numel(badChans)
                    ch = badChans(k);
                    lbl = safe_chan_label(EEG, ch);
                    logmsg(logf, '  - Ch %d (%s): %s', ch, lbl, reasons{k});
                end
            end

            % 4) Manual decision: which to interp (default = none)
            interpChans = prompt_channel_interp(EEG, badChans);

            if isempty(interpChans)
                logmsg(logf, '[INITREJ] No channels interpolated (namual decision).');
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
            else
                logmsg(logf, '[INITREJ] Interpolating channels (spherical): %s', vec2str(interpChans));
                EEG = pop_interp(EEG, interpChans, 'spherical');
                EEG = eeg_checkset(EEG);

<<<<<<< HEAD
=======
                % Suggested but not interpolated (and why suggested)
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
                notInterp = setdiff(badChans, interpChans);
                if ~isempty(notInterp)
                    logmsg(logf, '[INITREJ] Suggested but NOT interpolated: %s', vec2str(notInterp));
                    for k = 1:numel(notInterp)
                        ch = notInterp(k);
                        idx = find(badChans == ch, 1);
                        lbl = safe_chan_label(EEG, ch);
                        if ~isempty(idx)
                            logmsg(logf, '  - Ch %d (%s): %s', ch, lbl, reasons{idx});
<<<<<<< HEAD
=======
                        else
                            logmsg(logf, '  - Ch %d (%s): (reason missing)', ch, lbl);
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
                        end
                    end
                end
            end

            tags{end+1} = nextTag;
<<<<<<< HEAD
            save_stage(ST.INITREJ, P, subjid, tags, EEG, logf);
        end
    end

    %% ---------------
    % ICA (+ ICLabel)
    % ----------------
    if cfg.preproc.ica.enabled
        nextTag = char(string(cfg.preproc.ica.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.ICA, P, subjid, tags, nextTag, logf);
        if ~didLoad
            logmsg(logf, '[ICA] method = %s', char(string(cfg.preproc.ica.method)));

            [EEGtrain, segInfo] = make_ica_training_copy(EEG, cfg, logf);
            EEGtrain = pop_runica(EEGtrain, 'icatype', char(string(cfg.preproc.ica.method)));
            EEGtrain = eeg_checkset(EEGtrain);

            EEG.icaweights = EEGtrain.icaweights;
            EEG.icasphere = EEGtrain.icasphere;
            EEG.icawinv = EEGtrain.icawinv;
            EEG.icachansind = EEGtrain.icachansind;
            EEG = eeg_checkset(EEG, 'ica');

            trainedOn = 'FULL';
            if isfield(segInfo, 'removed') && logical(segInfo.removed)
                trainedOn = 'CLEAN-COPY';
            end

            nIntervals = NaN;
            pctTime = NaN;
            if isfield(segInfo, 'n_intervals'); nIntervals = segInfo.n_intervals; end
            if isfield(segInfo, 'pct_time'); pctTime = segInfo.pct_time; end
            
                logmsg(logf, '[ICA] Trained on %s. badseg_removed = %d intervals = %g pct = %.2f', trainedOn, logical(isfield(segInfo, 'removed') && segInfo.removed), nIntervals, pctTime);

            tags{end+1} = nextTag;

            % ICLabel suggest + QC + manual reject
=======
            save_stage(P.STAGE.INITREJ, P, subjid, tags, EEG, logf);
        end
    end


    % ------------------
    % ICA
    % ------------------
    if cfg.preproc.ica.enabled
        nextTag = char(string(cfg.preproc.ica.tag));
        [EEG, tags, didLoad] = maybe_load_stage(P.STAGE.ICA, P, subjid, tags, nextTag, logf);
        if ~didLoad
            logmsg(logf, '[ICA] method = %s', char(string(cfg.preproc.ica.method)));

            % ICA training on clean copy (optional bad segment removal)
            [EEGtrain, segInfo] = make_ica_training_copy(EEG, cfg, logf);

            EEGtrain = pop_runica(EEGtrain, 'icatype', char(string(cfg.preproc.ica.method)));
            EEGtrain = eeg_checkset(EEGtrain);

            % Copy ICA solution back to full EEG
            EEG.icaweights  = EEGtrain.icaweights;
            EEG.icasphere   = EEGtrain.icasphere;
            EEG.icawinv     = EEGtrain.icawinv;
            EEG.icachansind = EEGtrain.icachansind;
            EEG = eeg_checkset(EEG, 'ica');

            logmsg(logf, '[ICA] Trained on %s. badseg_removed = %d intervals = %d pct = %.2f', ...
                tern(segInfo.removed, 'CLEAN-COPY', 'FULL'), segInfo.removed, segInfo.n_intervals, segInfo.pct_time);

            tags{end+1} = nextTag;

            % -------------------------------------------
            % ICLabel: suggest + manual (no auto removal)
            % -------------------------------------------
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
            if isfield(cfg.preproc.ica, 'iclabel') && isfield(cfg.preproc.ica.iclabel, 'enabled') && cfg.preproc.ica.iclabel.enabled
                try
                    EEG = iclabel(EEG);
                    EEG = eeg_checkset(EEG);

                    thr = struct();
                    if isfield(cfg.preproc.ica.iclabel, 'thresholds')
                        thr = cfg.preproc.ica.iclabel.thresholds;
                    end
                    [suggestICs, icReasons] = iclabel_suggest_reject(EEG, thr);

                    if ~isfield(EEG, 'etc') || isempty(EEG.etc)
                        EEG.etc = struct();
                    end
                    EEG.etc.iclabel = struct();
                    EEG.etc.iclabel.suggestICs = suggestICs;
                    EEG.etc.iclabel.reasons = icReasons;

<<<<<<< HEAD
                    % Save QC packets -> QC folder
                    try
                        save_ic_qc_packets(QC, subjid, EEG, suggestICs);
                        logmsg(logf, '[ICQC] Saved IC QC packets to QC.');
=======
                    % Pre-mark suggestions (still no removal)
                    if ~isfield(EEG, 'reject') || isempty(EEG.reject)
                        EEG.reject = struct();
                    end
                    if ~isfield(EEG.reject, 'gcompreject') || isempty(EEG.reject.gcompreject)
                        EEG.reject.gcompreject = zeros(1, size(EEG.icaweights, 1));
                    end
                    EEG.reject.gcompreject(:) = 0;
                    if ~isempty(suggestICs)
                        EEG.reject.gcompreject(suggestICs) = 1;
                    end

                    if isempty(suggestICs)
                        logmsg(logf, '[ICLABEL] No ICs suggested for rejection.');
                    else
                        logmsg(logf, '[ICLABEL] Suggested ICs for rejection: %s', vec2str(suggestICs));
                        for k = 1:numel(icReasons)
                            logmsg(logf, '  - %s', icReasons{k});
                        end
                    end

                    % Save QC packets (topoplot + PSD + snippet) for suggested ICs
                    try
                        save_ic_qc_packets(P, subjid, EEG, suggestICs);
                        logmsg(logf, '[ICQC] Saved IC QC packets for suggested ICs.');
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
                    catch ME
                        logmsg(logf, '[WARN] IC QC packet generation failed: %s', ME.message);
                    end

<<<<<<< HEAD
                    % Optional tag for ICLabel pass (so filename captures
                    % it)
=======
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
                    if isfield(cfg.preproc.ica.iclabel, 'tag') && ~isempty(cfg.preproc.ica.iclabel.tag)
                        tags{end+1} = char(string(cfg.preproc.ica.iclabel.tag));
                    end

<<<<<<< HEAD
                    % Numeric mertrics -> LOGS
                    try
                        icMetrics = compute_ic_psd_metrics(EEG, suggestICs);
                        write_ic_metrics_csv(LOGS, subjid, icMetrics);
=======
                    % Numeric PSD metrics for suggested ICs (log + save)
                    try
                        icMetrics = compute_ic_psd_metrics(EEG, suggestICs);
                        write_ic_metrics_csv(P, subjid, icMetrics);
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
                        log_ic_metrics(logf, icMetrics);
                        logmsg(logf, '[ICMET] Logged PSD metrics + wrote CSV.');
                    catch ME
                        logmsg(logf, '[WARN] IC PSD metric computation failed: %s', ME.message);
                    end

<<<<<<< HEAD
                    % Manual reject (default none)
=======
                    % ---- Manual decision (default = remove NONE) ----
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
                    removedICs = prompt_ic_reject(suggestICs);

                    if isempty(removedICs)
                        logmsg(logf, '[ICREJ] No ICs removed.');
                    else
                        logmsg(logf, '[ICREJ] Removing ICs: %s', vec2str(removedICs));
<<<<<<< HEAD
                        EEG = pop_subcomp(EEG, removedICs, 0);
                        EEG = eeg_checkset(EEG);
                    end

=======
                        EEG = pop_subcomp(EEG, removedICs, 0); % 0 = don't plot
                        EEG = eeg_checkset(EEG);
                    end

                    % Log suggested vs removed vs kept
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
                    keptSuggested = setdiff(suggestICs, removedICs);
                    if ~isempty(keptSuggested)
                        logmsg(logf, '[ICREJ] Suggested but NOT removed: %s', vec2str(keptSuggested));
                        for k = 1:numel(keptSuggested)
                            ic = keptSuggested(k);
                            idx = find(suggestICs == ic, 1);
                            if ~isempty(idx)
                                logmsg(logf, '  - %s', icReasons{idx});
                            end
                        end
                    end

                catch ME
                    logmsg(logf, '[WARN] ICLabel failed: %s', ME.message);
                end
            else
                logmsg(logf, '[ICLABEL] disabled.');
            end
<<<<<<< HEAD

            save_stage(ST.ICA, P, subjid, tags, EEG, logf);
        end
    end

    %% -----
    % EPOCH
    % ------
    if cfg.preproc.epoch.enabled
        nextTag = char(string(cfg.preproc.epoch.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.EPOCH, P, subjid, tags, nextTag, logf);
        if ~didLoad
            ev = cfg.preproc.epoch.event_types;
            tmin = cfg.preproc.epoch.tmin_sec;
            tmax = cfg.preproc.epoch.tmax_sec;

            logmsg(logf, '[EPOCH] events = %s window = [%.3f %.3f] sec', strjoin(string(ev), ','), tmin, tmax);
            if isstring(ev); ev = cellstr(ev); end

            validate_events_before_epoch(EEG, ev, logf);

            EEG = pop_epoch(EEG, ev, [tmin tmax]);
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
            save_stage(ST.EPOCH, P, subjid, tags, EEG, logf);
        end
    end

    %% --------
    % BASELINE
    % ---------
    if cfg.preproc.baseline.enabled
        nextTag = char(string(cfg.preproc.baseline.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.BASE, P, subjid, tags, nextTag, logf);
        if ~didLoad
            win = cfg.preproc.baseline.window_sec;
            logmsg(logf, '[BASELINE] window = [%.3f %.3f] sec', win(1), win(2));

            EEG = pop_rmbase(EEG, win * 1000); % ms
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
            save_stage(ST.BASE, P, subjid, tags, EEG, logf);
        end
    end

    logmsg(logf, '===== PREPROC DONE sub-%03d =====', subjid);
end

=======

            save_stage(P.STAGE.ICA, P, subjid, tags, EEG, logf);
        end
    end

    % -----
    % EPOCH
    % -----
    if cfg.preproc.epoch.enabled
        nextTag = char(string(cfg.preproc.epoch.tag));
        [EEG, tags, didLoad] = maybe_load_stage(P.STAGE.EPOCH, P, subjid, tags, nextTag, logf);
        if ~didLoad
            ev = cfg.preproc.epoch.event_types;
            tmin = cfg.preproc.epoch.tmin_sec;
            tmax = cfg.preproc.epoch.tmax_sec;

            logmsg(logf, '[EPOCH] events = %s window = [%.3f %.3f] sec', strjoin(string(ev), ','), tmin, tmax);
            if isstring(ev); ev = cellstr(ev); end

            validate_events_before_epoch(EEG, ev, logf);

            EEG = pop_epoch(EEG, ev, [tmin tmax]);
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
            save_stage(P.STAGE.EPOCH, P, subjid, tags, EEG, logf);
        end
    end

    % --------
    % BASELINE
    % --------
    if cfg.preproc.baseline.enabled
        nextTag = char(string(cfg.preproc.baseline.tag));
        [EEG, tags, didLoad] = maybe_load_stage(P.STAGE.BASE, P, subjid, tags, nextTag, logf);
        if ~didLoad
            win = cfg.preproc.baseline.window_sec;
            logmsg(logf, '[BASELINE] window = [%.3f %.3f] sec', win(1), win(2));

            EEG = pop_rmbase(EEG, win * 1000); % ms
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
            save_stage(P.STAGE.BASE, P, subjid, tags, EEG, logf);
        end
    end

    logmsg(logf, '==== PREPROC DONE sub=%03d ====', subjid);
end

end

% ---------------------------
% HELPERS
% ---------------------------

function rawPath = resolve_raw_file(P, cfg, subjid)
rawPath = '';

pat = '';
if isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'pattern')
    pat = char(string(cfg.exp.raw.pattern));
end

if ~isempty(pat)
    fname = sprintf(pat, subjid);
    candidate = fullfile(P.INPUT.EXP, fname);
    if exist(candidate, 'file')
        rawPath = candidate;
        return;
    end
end

doRec = true;
if isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'search_recursive')
    doRec = logical(cfg.exp.raw.search_recursive);
end

if doRec
    pat1 = sprintf('*sub-%03d*.eeg', subjid);
    pat2 = sprintf('*%03d*.eeg', subjid);
    d = dir(fullfile(P.INPUT.EXP, '**', pat1));
    if isempty(d)
        d = dir(fullfile(P.INPUT.EXP, '**', pat2));
    end
    if ~isempty(d)
        rawPath = fullfile(d(1).folder, d(1).name);
    end
end
end

function save_stage(stageDir, P, subjid, tags, EEG, logf)
if ~exist(stageDir, 'dir'); mkdir(stageDir); end
fname = P.NAMING.fname(subjid, tags);
outPath = fullfile(stageDir, fname);
logmsg(logf, '  [SAVE] %s', outPath);

[fp, fn, ext] = fileparts(outPath);
pop_saveset(EEG, 'filename', [fn ext], 'filepath', fp); % #ok<NASGU>
end

% Montage application helper
function EEG = apply_montage_biosemi_from_csv(P, cfg, EEG, logf, subjid)
% Apply BioSemi A/B channel mapping using a CSV resource.
% - selects A1-32 and B1-32 if requested
% - relabels channels according to CSV
% - optionally applied .elp lookup for coords
%
% cfg.exp.montage.csv : filename under P.RESOURCE
% cfg.exp.montage.select_ab_only : true/false
% cfg.exp.montage.do_lookup : true/false

% --- Resolve CSV path ---
csvPath = "";
if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'csv')
    csvPath = string(cfg.exp.montage.csv);
end
if strlength(csvPath) == 0
    error('Montage enabled but cfg.exp.montage.csv is missing.');
end

csvPath2 = fullfile(P.RESOURCE, char(csvPath));
if isfile(csvPath2)
    csvPath = string(csvPath2);
end
if ~isfile(csvPath)
    error('Montage CSV not found: %s', char(csvPath));
end

% --- Optionally select only A/B scalp channels ---
selectAB = true;
if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'select_ab_only')
    selectAB = logical(cfg.exp.montage.select_ab_only);
end

if selectAB
    keep = [ ...
        arrayfun(@(x) sprintf('A%d', x), 1:32, 'UniformOutput', false), ...
        arrayfun(@(x) sprintf('B%d', x), 1:32, 'UniformOutput', false) ...
    ];
    EEG = pop_select(EEG, 'channel', keep);
    EEG = eeg_checkset(EEG);
    logmsg(logf, '[MONTAGE] Selected A1-32 & B1-32 only. nbchan = %d', EEG.nbchan);
end

% --- Read mapping ---
T = readtable(char(csvPath), 'Delimiter', ',', 'TextType', 'string');
T.Properties.VariableNames = lower(string(T.Properties.VariableNames));
reqCols = ["raw_label", "std_label"];
if ~all(ismember(reqCols, string(T.Properties.VariableNames)))
    error('Montage CSV must contain columns: raw_label, std_label');
end

raw = string(T.raw_label);
std = string(T.std_label);

% --- Validate mapping uniqueness ---
if numel(unique(raw)) ~= numel(raw)
    error('Montage CSV has duplicate raw_label entries.');
end
if numel(unique(std)) ~= numel(std)
    error('Montage CSV has duplicate std_label entries.');
end

% --- Apply relabelling by matching existing labels ---
cur = string({EEG.chanlocs.labels});
curU = upper(cur);
rawU = upper(raw);

missing = setdiff(rawU, curU);
if ~isempty(missing)
    error('Montage CSV expects channels not present in EEG: %s', strjoin(missing, ', '));
end

for i = 1:numel(rawU)
    idx = find(curU == rawU(i), 1);
    EEG.chanlocs(idx).labels = char(std(i));
end
EEG = eeg_checkset(EEG);

% --- Final duplicate check ---
labs = string({EEG.chanlocs.labels});
if numel(unique(upper(labs))) ~= numel(labs)
    error('Relabel produced duplicate labels (case-insensitive). Check montage CSV.');
end

% --- Sanity-check expected midline ---
mustHave = upper(["FPZ","AFZ","FZ","FCZ","CZ","CPZ","PZ","POZ","OZ","IZ"]);
have = upper(labs);
missingStd = setdiff(mustHave, have);
if ~isempty(missingStd)
    logmsg(logf, '[WARN] Montage missing expected midline labels: %s', strjoin(missingStd, ', '));
end

logmsg(logf, '[MONTAGE] Relabel complete from CSV: %s', char(csvPath));

% --- Optional ELP lookup ---
doLookup = true;
if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'do_lookup')
    doLookup = logical(cfg.exp.montage.do_lookup);
end

if doLookup
    if isfield(P, 'CORE') && isfield(P.CORE, 'ELP_FILE') && exist(P.CORE.ELP_FILE, 'file')
        EEG = pop_chanedit(EEG, 'lookup', P.CORE.ELP_FILE);
        EEG = eeg_checkset(EEG);
        logmsg(logf, '[MONTAGE] Applied coord lookup from ELP: %s', P.CORE.ELP_FILE);
    else
        logmsg(logf, '[WARN] ELP file missing; skipping lookup (P.CORE.ELP_FILE).');
    end
end

write_channelmap_tsv(P, subjid, EEG, table(raw, std));
end

% Channel map audit writer
function write_channelmap_tsv(P, subjid, EEG)
% Writes an audit file of channel relabeling:
% columns: index, final_label

outPath = fullfile(P.STAGE.LOGS, sprintf('sub-%03d_channelmap_applied.tsv', subjid));
fid = fopen(outPath, 'w');
if fid < 0
    warning('Could not write channel map TSV: %s', outPath);
    return;
end

fprintf(fid, "index\tlabel\n");
for i = 1:EEG.nbchan
    lbl = EEG.chanlocs(i).labels;
    fprintf(fid, "%d\t%s\n", i, lbl);
end

fclose(fid);
end

function logmsg(fid, fmt, varargin)
ts = datetime("now", "Format", "yyyy-MM-dd HH:mm:ss");
msg = sprintf(fmt, varargin{:});
fprintf(fid, '[%s] %s\n', char(ts), msg);
if fid ~= 1
    fprintf(1, '[%s] %s\n', char(ts), msg);
end
end

function safeClose(fid)
if fid ~= 1 && fid > 0
    fclose(fid);
end
end

% Resume/Skip helper
function [EEG, tags, didLoad] = maybe_load_stage(stageDir, P, subjid, tags, nextTag, logf)
didLoad = false;
tags2 = tags;
if ~isempty(nextTag)
    tags2{end+1} = nextTag;
end
fname = P.NAMING.fname(subjid, tags2);
fpath = fullfile(stageDir, fname);
if exist(fpath, 'file')
    logmsg(logf, '[SKIP] Found existing stage file, loading: %s', fpath);
    EEG = pop_loadset('filename', fname, 'filepath', stageDir);
    EEG = eeg_checkset(EEG);
    tags = tags2;
    didLoad = true;
end
end

% Channel label normalization
function EEG = normalize_chan_labels(EEG)
if ~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs)
    return;
end
for i = 1:numel(EEG.chanlocs)
    if isfield(EEG.chanlocs(i), 'labels') && ~isempty(EEG.chanlocs(i).labels)
        lbl = string(EEG.chanlocs(i).labels);
        lbl = strtrim(lbl);
        lbl = replace(lbl, " ", "");
        lbl = replace(lbl, "-", "");
        EEG.chanlocs(i).labels = char(lbl);
    end
end
end

% ---- INITREJ: suggestions + plots + manual prompt ----

function [badChans, reasons, metrics] = suggest_bad_channels(EEG)
% Conservative automated suggestions using EEGLAB pop_rejchan (prob +
% kurtosis)
% Also returns metrics for plotting (STD/RMS)

badChans = [];
reasons = {};

metrics = struct();
metrics.chan_rms = sqrt(mean(double(EEG.data).^2, 2)); % [nChan x 1]
metrics.chan_std = std(double(EEG.data), 0, 2);

% Probability-based suggestions
try
    [~, badP] = pop_rejchan(EEG, 'threshold', 5, 'norm', 'on', 'measure', 'prob');
    if ~isempty(badP)
        for c = badP(:)'
            badChans(end+1) = c; %#ok<AGROW>
            reasons{end+1} = 'probability z>5 (pop_rejchan)'; %#ok<AGROW>
        end
    end
catch
end

% Kurtosis-based suggestions
try
    [~, badK] = pop_rejchan(EEG, 'threshold', 5, 'norm', 'on', 'measure', 'kurt');
    if ~isempty(badK)
        for c = badK(:)'
            if ~ismember(c, badChans)
                badChans(end+1) = c; %#ok<AGROW>
                reasons{end+1} = 'kurtosis z>5 (pop_rejchan)'; %#ok<AGROW>
            else
                idx = find(badChans == c, 1);
                reasons{idx} = [reasons{idx} ' + kurtosis z>5'];
            end
        end
    end
catch
end

[badChans, sortIdx] = sort(badChans);
reasons = reasons(sortIdx);
end

function make_initrej_plots(P, subjid, EEG, metrics, badChans)
% Save diagnostic plots to logs/:
%   - histogram of channel std/rms
%   - barplot of channel std/rms with suggested chans marked
%   - topoplot of std/rms (if chanlocs present)
%   - channel PSD overview (median + IQR)
%   - PSD overlay for suggested channels

outDir = P.STAGE.LOGS;

save_hist(outDir, subjid, metrics.chan_std, 'chan_std');
save_hist(outDir, subjid, metrics.chan_rms, 'chan_rms');

save_bar(outDir, subjid, metrics.chan_std, badChans, 'chan_std');
save_bar(outDir, subjid, metrics.chan_rms, badChans, 'chan_rms');

if has_chanlocs(EEG)
    save_topo_metric(outDir, subjid, EEG, metrics.chan_std, 'STD');
    save_topo_metric(outDir, subjid, EEG, metrics.chan_rms, 'RMS');
end

save_channel_psd_overview(outDir, subjid, EEG);

if ~isempty(badChans)
    save_channel_psd_badchans(outDir, subjid, EEG, badChans);
end

% write label index map for quick referencing
fid = fopen(fullfile(outDir, sprintf('sub-%03d_chanlabels.txt', subjid)), 'w');
if fid > 0
    labs = get_chan_labels(EEG);
    for k = 1:numel(labs)
        fprintf(fid, '%d\t%s\n', k, labs{k});
    end
    fclose(fid);
end
end

function save_hist(outDir, subjid, v, name)
h = figure('Visible', 'off');
histogram(v);
title(sprintf('sub=%03d %s histogram', subjid, name), 'Interpreter', 'none');
xlabel(name); ylabel('Count');
saveas(h, fullfile(outDir, sprintf('sub-%03d_initrej_hist_%s.png', subjid, name)));
close(h);
end

function save_bar(outDir, subjid, v, badChans, name)
h = figure('Visible', 'off');
bar(v);
title(sprintf('sub-%03d %s (suggested marked)', subjid, name), 'Interpreter', 'none');
xlabel('Channel index'); ylabel(name);
hold on;
for c = badChans(:)'
    xline(c, '--');
end
hold off;
saveas(h, fullfile(outDir, sprintf('sub-%03d_initrej_bar_%s.png', subjid, name)));
close(h);
end

function tf = has_chanlocs(EEG)
tf = isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs(1), 'X');
end

function save_topo_metric(outDir, subjid, EEG, v, label)
h = figure('Visible', 'off');
topoplot(v, EEG.chanlocs, 'electrodes', 'on');
title(sprintf('sub-%03d channel %s topoplot', subjid, label), 'Interpreter', 'none');
saveas(h, fullfile(outDir, sprintf('sub-%03d_initrej_topo_%s.png', subjid, lower(label))));
close(h);
end

function save_channel_psd_overview(outDir, subjid, EEG)
fs = EEG.srate;
data = double(EEG.data);
nChan = size(data, 1);

win = round(fs*2);
nover = round(win*0.5);
nfft = max(2^nextpow2(win), win);

Pxx = zeros(nfft/2+1, nChan);
for ch = 1:nChan
    [pxx, f] = pwelch(data(ch,:), win, nover, nfft, fs);
    Pxx(:, ch) = 10*log10(pxx);
end

med = median(Pxx, 2);
q1 = prctile(Pxx, 25, 2);
q3 = prctile(Pxx, 75, 2);

h = figure('Visible', 'off');
plot(f, med); hold on;
plot(f, q1, ':'); plot(f, q3, ':');
xlim([0 80]);
xlabel('Hz'); ylabel('Power (dB)');
title(sprintf('sub-%03d channel PSD (median + IQR)', subjid), 'Interpreter', 'none');
legend({'Median', '25%', '75%'}, 'Location', 'northeast');
saveas(h, fullfile(outDir, sprintf('sub-%03d_initrej_psd_overview.png', subjid)));
close(h);
end

function save_channel_psd_badchans(outDir, subjid, EEG, badChans)
fs = EEG.srate;
data = double(EEG.data);

win = round(fs*2);
nover = round(win*0.5);
nfft = max(2^nextpow2(win), win);

h = figure('Visible', 'off'); hold on;
for ch = badChans(:)'
    [pxx, f] = pwelch(data(ch, :), win, nover, nfft, fs);
    plot(f, 10*log10(pxx));
end
xlim([0 80]);
xlabel('Hz'); ylabel('Power (dB)');
title(sprintf('sub-%03d PSD overlay: suggested channels %s', subjid, vec2str(badChans)), 'Interpreter', 'none');
saveas(h, fullfile(outDir, sprintf('sub-%03d_initrej_psd_badchans.png', subjid)));
close(h);
end

function interpChans = prompt_channel_interp(EEG, suggested)
fprintf('\n[INITREJ] Suggested bad channels: %s\n', vec2str(suggested));
if ~isempty(suggested)
    for k = 1:numel(suggested)
        ch = suggested(k);
        fprintf('   %d) %d (%s)\n', k, ch, safe_chan_label(EEG, ch));
    end
end
fprintf('\nType channel indices to INTERPOLATE (e.g., [1 2 17])\n');
fprintf('Default = none (press Enter or type []).\n');
resp = input('Channels to interpolate: ', 's');
resp = strtrim(resp);

if isempty(resp) || strcmp(resp, '[]')
    interpChans = [];
    return;
end

interpChans = str2num(resp); %#ok<ST2NM>
if ~isnumeric(interpChans)
    error('Manual interpolation must be numeric vector like [1 2 3] or [].');
end
interpChans = unique(interpChans(:))';
end

function lbl = safe_chan_label(EEG, ch)
lbl = '';
try
    if isfield(EEG, 'chanlocs') && numel(EEG.chanlocs) >= ch && isfield(EEG.chanlocs(ch), 'labels')
        lbl = EEG.chanlocs(ch).labels;
    end
catch
end
if isempty(lbl)
    lbl = sprintf('Ch%d', ch);
end
end

function labels = get_chan_labels(EEG)
n = EEG.nbchan;
labels = cell(n, 1);
for i = 1:n
    labels{i} = safe_chan_label(EEG, i);
end
end

function s = vec2str(v)
if isempty(v), s = '[]'; return; end
s = ['[' sprintf('%d ', v) ']'];
s = strrep(s, ' ]', ']');
end

% Per-channel PSD metrics + CSV + topo
function chanPSD = compute_channel_psd_metrics(EEG)
fs = EEG.srate;
data = double(EEG.data);
nChan = size(data, 1);

win = round(fs*2);
nover = round(win*0.5);
nfft = max(2^nextpow2(win), win);

chanPSD = struct();
chanPSD.line_ratio  = zeros(nChan, 1);
chanPSD.hf_ratio    = zeros(nChan, 1);
chanPSD.drift_ratio = zeros(nChan, 1);
chanPSD.alpha_ratio = zeros(nChan, 1);

for ch = 1:nChan
    [pxx, f] = pwelch(data(ch,:), win, nover, nfft, fs);

    line = bp_psd(f, pxx, [59 61]);
    lo   = bp_psd(f, pxx, [55 59]);
    hi   = bp_psd(f, pxx, [61 65]);
    chanPSD.line_ratio(ch) = line / max(lo + hi, eps);

    hf = bp_psd(f, pxx, [20 40]);
    lf = bp_psd(f, pxx, [1 12]);
    chanPSD.hf_ratio(ch) = hf / max(lf, eps);

    drift = bp_psd(f, pxx, [1 2]);
    chanPSD.drift_ratio(ch) = drift / max(lf, eps);

    a = bp_psd(f, pxx, [8 12]);
    allp = bp_psd(f, pxx, [1 40]);
    chanPSD.alpha_ratio(ch) = a / max(allp, eps);
end
end

function p = bp_psd(f, pxx, band)
m = (f >= band(1) & f < band(2));
if ~any(m); p = 0; return; end
p = trapz(f(m), pxx(m));
end

function write_channel_psd_csv(P, subjid, EEG, chanPSD)
outDir = P.STAGE.LOGS;
csvPath = fullfile(outDir, sprintf('sub-%03d_chan_psd_metrics.csv', subjid));
fid = fopen(csvPath, 'w');
if fid < 0; return; end
fprintf(fid, 'chan_idx,label,line_ratio,hf_ratio,drift_ratio,alpha_ratio\n');
for ch = 1:EEG.nbchan
    fprintf(fid, '%d,%s,%.6f,%.6f,%.6f,%.6f\n', ch, safe_chan_label(EEG,ch), ...
        chanPSD.line_ratio(ch), chanPSD.hf_ratio(ch), chanPSD.drift_ratio(ch), chanPSD.alpha_ratio(ch));
end
fclose(fid);
end

function save_chan_psd_topos(P, subjid, EEG, chanPSD)
if ~has_chanlocs(EEG); return; end
outDir = P.STAGE.LOGS;
save_topo_metric(outDir, subjid, EEG, chanPSD.line_ratio, 'LINE_RATIO');
save_topo_metric(outDir, subjid, EEG, chanPSD.hf_ratio, 'HF_RATIO');
save_topo_metric(outDir, subjid, EEG, chanPSD.drift_ratio, 'DRIFT_RATIO');
save_topo_metric(outDir, subjid, EEG, chanPSD.alpha_ratio, 'ALPHA_RATIO');
end

% ICA training copy with optional bad segment removal
function [EEGtrain, segInfo] = make_ica_training_copy(EEG, cfg, logf)
EEGtrain = EEG;
segInfo = struct('removed', false, 'n_intervals', 0, 'pct_time', 0, 'intervals', []);
if ~isfield(cfg, 'preproc') || ~isfield(cfg.preproc, 'initrej') || ~isfield(cfg.preproc.initrej, 'badseg') || ~cfg.preproc.initrej.badseg.enabled
    logmsg(logf, '[ICA-TRAIN] badseg disabled; using full data for ICA.');
    return;
end
thr = cfg.preproc.initrej.badseg.threshold_uv;
logmsg(logf, '[ICA-TRAIN] Detecting bad segments (thr = %.1f uV) for ICA training copy.', thr);
x = double(EEG.data);
badSamp = any(abs(x) > thr, 1);
intervals = mask_to_intervals(badSamp);
segInfo.intervals = intervals;
segInfo.n_intervals = size(intervals, 1);
segInfo.pct_time = 100 * (sum(badSamp) / numel(badSamp));
if isempty(intervals)
    logmsg(logf, '[ICA-TRAIN] No bad segments detected; using full data.');
    return;
end
logmsg(logf, '[ICA-TRAIN] Detected %d intervals (%.2f%% of samples).', segInfo.n_intervals, segInfo.pct_time);
doRemove = prompt_yesno('Remove detected bad segments from ICA training copy? (y/n) [n]: ', false);
if ~doRemove
    logmsg(logf, '[ICA-TRAIN] Keeping all segments for ICA training (manual decision).');
    return;
end
EEGtrain = pop_select(EEGtrain, 'nopoint', intervals);
EEGtrain = eeg_checkset(EEGtrain);
segInfo.removed = true;
logmsg(logf, '[ICA-TRAIN] Removed bad segments from ICA training copy.');
end

function intervals = mask_to_intervals(mask)
mask = mask(:)';
if ~any(mask)
    intervals = [];
    return;
end
d = diff([0 mask 0]);
starts = find(d == 1);
ends   = find(d == -1) - 1;
intervals = [starts(:) ends(:)];
end

function tf = prompt_yesno(prompt, defaultTF)
resp = input(prompt, 's');
resp = lower(strtrim(resp));
if isempty(resp)
    tf = defaultTF;
elseif any(strcmp(resp, {'y', 'yes'}))
    tf = true;
elseif any(strcmp(resp, {'n', 'no'}))
    tf = false;
else
    tf = defaultTF;
end
end

function s = tern(cond, a, b)
if cond
    s = a;
else
    s = b;
end
end

% Event validation before epoching
function validate_events_before_epoch(EEG, wanted, logf)
if isstring(wanted)
    wanted = cellstr(wanted);
end
if isempty(EEG.event)
    error('No EEG.event present; cannot epoch.');
end
types = {EEG.event.type};
typesStr = strings(1, numel(types));
for k = 1:numel(types)
    typesStr(k) = string(types{k});
end
u = unique(typesStr);
logmsg(logf, '[EVENTS] Unique event types present (%d): %s', numel(u), strjoin(u, ', '));
for k = 1:numel(u)
    c = sum(typesStr == u(k));
    logmsg(logf, '  [EVENTS] %s: %d', u(k), c);
end
wantedStr = string(wanted);
hit = intersect(u, wantedStr);
if isempty(hit)
    error('None of requested event_types found: %s', strjoin(wantedStr, ', '));
else
    for k = 1:numel(hit)
        logmsg(logf, '[EVENTS] Will epoch on "%s" (n = %d)', hit(k), sum(typesStr == hit(k)));
    end
end
end

% ---- ICLabel: suggest + QC packets + manual rejection ----
function [suggestICs, reasons] = iclabel_suggest_reject(EEG, thr)
% Uses EEG.etc.ic_classification.ICLabel.classifications (Nx7)
% Order: [Brain Muscle Eye Heart LineNoise ChannelNoise Other]

suggestICs = [];
reasons = {};

if ~isfield(EEG, 'etc') || ~isfield(EEG.etc, 'ic_classification') || ...
        ~isfield(EEG.etc.ic_classification, 'ICLabel') || ...
        ~isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
    return;
end

C = EEG.etc.ic_classification.ICLabel.classifications;
if isempty(C)
    return;
end

for ic = 1:size(C, 1)
    pBrain  = C(ic, 1);
    pMus    = C(ic, 2);
    pEye    = C(ic, 3);
    pHeart  = C(ic, 4);
    pLine   = C(ic, 5);
    pChan   = C(ic, 6);
    pOther  = C(ic, 7);

    hits = {};
    if isfield(thr,'eye') && pEye >= thr.eye
        hits{end+1} = sprintf('eye=%.2f>=%.2f', pEye, thr.eye); %#ok<AGROW>
    end
    if isfield(thr,'muscle') && pMus >= thr.muscle
        hits{end+1} = sprintf('muscle=%.2f>=%.2f', pMus, thr.muscle); %#ok<AGROW>
    end
    if isfield(thr,'heart') && pHeart >= thr.heart
        hits{end+1} = sprintf('heart=%.2f>=%.2f', pHeart, thr.heart); %#ok<AGROW>
    end
    if isfield(thr,'line_noise') && pLine >= thr.line_noise
        hits{end+1} = sprintf('line=%.2f>=%.2f', pLine, thr.line_noise); %#ok<AGROW>
    end
    if isfield(thr,'channel_noise') && pChan >= thr.channel_noise
        hits{end+1} = sprintf('chanNoise=%.2f>=%.2f', pChan, thr.channel_noise); %#ok<AGROW>
    end

    if ~isempty(hits)
        suggestICs(end+1) = ic; %#ok<AGROW>
        reasons{end+1} = sprintf('IC %d: %s (brain = %.2f other = %.2f)', ...
            ic, strjoin(hits, ', '), pBrain, pOther); %#ok<AGROW>
    end
end
end

function save_ic_qc_packets(P, subjid, EEG, icList)
% Saves per-IC QC figures for suggested ICs:
%   - scalp map (icawinv)
%   - IC PSD (Welch)
%   - IC activation snippet

if isempty(icList); return; end

outDir = fullfile(P.STAGE.LOGS, sprintf('sub-%03d_icqc', subjid));
if ~exist(outDir, 'dir'); mkdir(outDir); end

fs = EEG.srate;

W = EEG.icaweights * EEG.icasphere;
X = double(EEG.data(EEG.icachansind, :));
act = W * X; % [nIC x time]

C = [];
if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification') && ...
        isfield(EEG.etc.ic_classification, 'ICLabel') && ...
        isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
    C = EEG.etc.ic_classification.ICLabel.classifications;
end

for ic = icList(:)'
    h = figure('Visible', 'off');

    % 1) Scalp map
    subplot(2,2,1);
    if isfield(EEG, 'icawinv') && has_chanlocs(EEG)
        topoplot(EEG.icawinv(:,ic), EEG.chanlocs, 'electrodes', 'on');
        title(sprintf('IC %d scalp map', ic), 'Interpreter', 'none');
    else
        axis off; title('No chanlocs/icawinv');
    end

    % 2) PSD
    subplot(2,2,2);
    win = round(fs * 2);
    nover = round(win * 0.5);
    nfft = max(2^nextpow2(win), win);
    [pxx, f] = pwelch(act(ic,:), win, nover, nfft, fs);
    plot(f, 10*log10(pxx));
    xlim([0 80]);
    xlabel('Hz'); ylabel('Power (dB)');
    title('IC PSD', 'Interpreter', 'none');

    % 3/4) time series snippet
    subplot(2, 2, [3 4]);
    nSamp = min(size(act, 2), round(fs * 10));
    plot((0:nSamp-1)/fs, act(ic,1:nSamp));
    xlabel('Time (s)'); ylabel('a.u.');
    title('IC activation (first 10s)', 'Interpreter', 'none');

    % Title with ICLabel probs if available
    if ~isempty(C) && size(C, 1) >= ic
        p = C(ic,:);
        ttl = sprintf('sub-%03d IC %d | B%.2f M%.2f E%.2f H%.2f L%.2f C%.2f O%.2f', ...
            subjid, ic, p(1), p(2), p(3), p(4), p(5), p(6), p(7));
    else
        ttl = sprintf('sub-%03d IC %d', subjid, ic);
    end
    sgtitle(ttl, 'Interpreter', 'none');

    saveas(h, fullfile(outDir, sprintf('sub-%03d_ic%03d_qc.png', subjid, ic)));
    close(h);
end
end

function removedICs = prompt_ic_reject(suggestICs)
fprintf('\n[ICREJ] ICLabel suggested ICs: %s\n', vec2str(suggestICs));
fprintf('Review QC figs in logs/sub-XXX_icqc/ before deciding.\n');
fprintf('Type IC indices to REMOVE (e.g., [1 3 7]).\n');
fprintf('Default = remove none (press Enter or type []).\n');

resp = input('ICs to remove: ', 's');
resp = strtrim(resp);

if isempty(resp) || strcmp(resp, '[]')
    removedICs = [];
    return;
end

removedICs = str2num(resp); %#ok<ST2NM>
if isempty(removedICs)
    removedICs = [];
else
    removedICs = unique(removedICs(:))';
end
end

function icMetrics = compute_ic_psd_metrics(EEG, icList)
% Returns struct array with per-IC PSD summary metrics

if isempty(icList)
    icMetrics = struct('ic', {}, 'peak_hz', {}, 'bp', {}, 'hf_ratio', {}, 'line_ratio', {});
    return;
end

fs = EEG.srate;

% activations
W = EEG.icaweights * EEG.icasphere;
X = double(EEG.data(EEG.icachansind, :));
act = W * X; % [nIC x time]

win = round(fs * 2);
nover = round(win * 0.5);
nfft = max(2^nextpow2(win), win);

% Bands (Hz)
bands = struct( ...
    'delta', [1 4], ...
    'theta', [4 8], ...
    'alpha', [8 12], ...
    'beta',  [13 30], ...
    'gamma', [30 45]);

icMetrics = repmat(struct('ic',NaN,'peak_hz',NaN,'bp',struct(),'hf_ratio',NaN,'line_ratio',NaN), numel(icList), 1);

for k = 1:numel(icList)
    ic = icList(k);

    [pxx, f] = pwelch(act(ic,:), win, nover, nfft, fs); % linear power

    % Peak frequency (0.5-40 Hz)
    bandMask = (f >= 0.5 & f <= 40);
    pband = pxx(bandMask);
    fband = f(bandMask);
    [~, idx] = max(pband);
    peakHz = fband(idx);

    % Bandpower (integral)
    bp = struct();
    fn = fieldnames(bands);
    for j = 1:numel(fn)
        b = bands.(fn{j});
        m = (f >= b(1) & f < b(2));
        bp.(fn{j}) = trapz(f(m), pxx(m));
    end

    % HF ratio: (20-40)/(1-12)
    hf = bandpower_from_psd(f, pxx, [20 40]);
    lf = bandpower_from_psd(f, pxx, [1 12]);
    hf_ratio = hf / max(lf, eps);

    % Line noise ratio: (59-61)/(55-65 excluding 59-61)
    line = bandpower_from_psd(f, pxx, [59 61]);
    lo   = bandpower_from_psd(f, pxx, [55 59]);
    hi   = bandpower_from_psd(f, pxx, [61 65]);
    denom = max(lo + hi, eps);
    line_ratio = line / denom;

    icMetrics(k).ic = ic;
    icMetrics(k).peak_hz = peakHz;
    icMetrics(k).bp = bp;
    icMetrics(k).hf_ratio = hf_ratio;
    icMetrics(k).line_ratio = line_ratio;
end
end

function p = bandpower_from_psd(f, pxx, band)
m = (f >= band(1) & f < band(2));
if ~any(m), p = 0; return; end
p = trapz(f(m), pxx(m));
end

function write_ic_metrics_csv(P, subjid, icMetrics)
outDir = P.STAGE.LOGS;
if ~exist(outDir, 'dir'); mkdir(outDir); end

csvPath = fullfile(outDir, sprintf('sub-%03d_ic_psd_metrics.csv', subjid));
fid = fopen(csvPath, 'w');
if fid < 0
    warning('Could not write IC metrics CSV: %s', csvPath);
    return;
end

% header
fprintf(fid, 'ic,peak_hz,delta,theta,alpha,beta,gamma,hf_ratio,line_ratio\n');

for k = 1:numel(icMetrics)
    bp = icMetrics(k).bp;
    fprintf(fid, '%d,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6f,%.6f\n', ...
        icMetrics(k).ic, icMetrics(k).peak_hz, ...
        bp.delta, bp.theta, bp.alpha, bp.beta, bp.gamma, ...
        icMetrics(k).hf_ratio, icMetrics(k).line_ratio);
end

fclose(fid);
end

function log_ic_metrics(logf, icMetrics)
if isempty(icMetrics); return; end
logmsg(logf, '[ICMET] Per-IC PSD summary (suggested ICs):');
for k = 1:numel(icMetrics)
    bp = icMetrics(k).bp;
    logmsg(logf, '  IC %d | peak = %.2f Hz | d = %.2e t = %.2e a = %.2e b = %.2e g = %.2e | HP = %.3f | line = %.3f', ...
        icMetrics(k).ic, icMetrics(k).peak_hz, ...
        bp.delta, bp.theta, bp.alpha, bp.beta, bp.gamma, ...
        icMetrics(k).hf_ratio, icMetrics(k).line_ratio);
end
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
end