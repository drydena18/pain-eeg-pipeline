function preproc_core(P, cfg)
% PREPROC_CORE Execute preprocessing according to cfg (from JSON)
% V 2.1.0
%
% V2.1.0 adds Stage 00 — multi-session EEG concatenation.
%   When cfg.preproc.concat.enabled = true, all session raw files for the
%   subject are loaded, merged with pop_mergeset (which adjusts event
%   latencies automatically across sessions), and saved to 00_concat/.
%   All downstream stages (01-08, 09_hilbert) are completely unaffected.
%   Single-session experiments (concat.enabled = false, the default) follow
%   the original code path with zero behavioural change.
%
%   Relevant cfg fields (safe defaults in normalize_preproc_defaults):
%     cfg.preproc.concat.enabled          true | false  (default false)
%     cfg.preproc.concat.n_sessions       integer       (default 1)
%     cfg.preproc.concat.session_pattern  printf string relative to P.INPUT.EXP
%                                         args: (subjid, sess) or (subjid, sess, subjid, sess)
%                                         Leave blank to use BIDS ses-NN layout or
%                                         the recursive alphabetical fallback.
%
%   Note on trial ordering: pop_mergeset appends sessions in order 1..N.
%   The participants_singletrial CSV (loaded by spec_load_singletrial_ratings)
%   must be sorted by (session, trial) for ratings to align correctly.
%
% Folder scheme (per subject):
%   PROJ_ROOT/<exp_out>/sub-XXX/00_concat  (multi-session only)
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
    ST.CONCAT   = fullfile(subRoot, '00_concat');
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

    % ================================================================
    % STAGE 00 — Multi-session concatenation  (or single-session load)
    %
    % After this block EEG is a continuous, montage-applied,
    % channel-location-loaded dataset ready for stage 01.
    % tags is {} (single-session) or {'concat'} (multi-session).
    % ================================================================
    doConcat = isfield(cfg.preproc, 'concat') && ...
               isfield(cfg.preproc.concat, 'enabled') && ...
               logical(cfg.preproc.concat.enabled);

    if doConcat
        % ----------------------------------------------------------
        % Multi-session path
        % ----------------------------------------------------------
        nSess     = cfg.preproc.concat.n_sessions;
        concatTag = 'concat';
        logmsg(logf, '[CONCAT] Multi-session mode: %d sessions.', nSess);

        [EEG, tags, didLoad] = maybe_load_stage(ST.CONCAT, P, subjid, tags, concatTag, logf);

        if ~didLoad
            try
                sessPaths = resolve_raw_session_files(P, cfg, subjid, nSess, logf);
            catch ME
                logmsg(logf, '[CONCAT][ERROR] %s -- skipping subject.', ME.message);
                continue;
            end

            sessEEGs = cell(nSess, 1);
            for s = 1:nSess
                logmsg(logf, '[CONCAT] Loading session %d: %s', s, sessPaths{s});
                Etmp = pop_biosig(sessPaths{s});
                Etmp = eeg_checkset(Etmp);
                Etmp = normalize_chan_labels(Etmp);

                % Montage applied per-session so channel counts match before merge
                if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'enabled') && ...
                        cfg.exp.montage.enabled
                    Etmp = apply_montage_biosemi_from_csv(P, cfg, Etmp, logf, subjid, LOGS);
                    Etmp = eeg_checkset(Etmp);
                end

                % Channel locations per-session (ensures chanlocs present on all)
                if isfield(cfg.exp, 'channel_locs') && isfield(cfg.exp.channel_locs, 'use_elp') && ...
                        cfg.exp.channel_locs.use_elp
                    elpPath = "";
                    if isfield(P, 'CORE') && isfield(P.CORE, 'ELP_FILE')
                        elpPath = string(P.CORE.ELP_FILE);
                    end
                    if strlength(elpPath) > 0 && exist(elpPath, 'file')
                        try
                            Etmp = pop_chanedit(Etmp, 'lookup', char(elpPath));
                            Etmp = eeg_checkset(Etmp);
                        catch ME
                            logmsg(logf, '[CONCAT][WARN] Chanlocs load failed sess %d: %s', s, ME.message);
                        end
                    end
                end

                logmsg(logf, '[CONCAT] Session %d: %d chans  %.1f s  srate=%.0f Hz', ...
                    s, Etmp.nbchan, Etmp.xmax, Etmp.srate);
                sessEEGs{s} = Etmp;
            end

            % Merge: pop_mergeset adjusts event latencies automatically
            logmsg(logf, '[CONCAT] Merging %d sessions with pop_mergeset...', nSess);
            EEG = sessEEGs{1};
            for s = 2:nSess
                EEG = pop_mergeset(EEG, sessEEGs{s});
            end
            EEG = eeg_checkset(EEG);

            logmsg(logf, '[CONCAT] Merged: %d chans  %.1f s  %d events', ...
                EEG.nbchan, EEG.xmax, length(EEG.event));

            % Provenance record
            EEG.etc.concat = struct();
            EEG.etc.concat.n_sessions    = nSess;
            EEG.etc.concat.session_files = sessPaths;

            tags{end+1} = concatTag;
            save_stage(ST.CONCAT, P, subjid, tags, EEG, logf);
        end

    else
        % ----------------------------------------------------------
        % Single-session path  (original behaviour, unchanged)
        % ----------------------------------------------------------
        rawPath = resolve_raw_file(P, cfg, subjid);
        if isempty(rawPath) || ~exist(rawPath, 'file')
            logmsg(logf, '[WARN] Raw file not found for sub-%03d. Skipping.', subjid);
            continue;
        end
        logmsg(logf, 'Raw File: %s', rawPath);

        rawPath = char(string(rawPath));
        assert(exist(rawPath, 'file') == 2, 'Raw file missing: %s', rawPath);

        logmsg(logf, 'RawPath class = %s | isstring = %d', class(rawPath), isstring(rawPath));

        EEG = pop_biosig(rawPath);
        EEG = eeg_checkset(EEG);
        EEG = normalize_chan_labels(EEG);

        % Montage (optional)
        if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'enabled') && cfg.exp.montage.enabled
            EEG = apply_montage_biosemi_from_csv(P, cfg, EEG, logf, subjid, LOGS);
            EEG = eeg_checkset(EEG);
        end

        % Channel Locations (.elp) if requested
        if isfield(cfg.exp, 'channel_locs') && isfield(cfg.exp.channel_locs, 'use_elp') && ...
                cfg.exp.channel_locs.use_elp
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
        end
    end

    % Deterministic seed per subject (reproducible ICA)
    rng(double(subjid), 'twister');
    logmsg(logf, '[RNG] rng(%d, twister)', subjid);

    %% ------
    % FILTER
    % -------
    if cfg.preproc.filter.enabled
        nextTag = char(string(cfg.preproc.filter.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.FILTER, P, subjid, tags, nextTag, logf, EEG);
        if ~didLoad
            hp = cfg.preproc.filter.highpass_hz;
            lp = cfg.preproc.filter.lowpass_hz;
            logmsg(logf, '[FILTER] %s hp = %.3f lp = %.3f', string(cfg.preproc.filter.type), hp, lp);

            EEG = pop_eegfiltnew(EEG, hp, lp);
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
            save_stage(ST.FILTER, P, subjid, tags, EEG, logf);
        end
    end

    %% ----------------------------
    % NOTCH (bandstop via revfilt)
    % -----------------------------
    if cfg.preproc.notch.enabled
        nextTag = char(string(cfg.preproc.notch.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.NOTCH, P, subjid, tags, nextTag, logf, EEG);
        if ~didLoad
            f0 = cfg.preproc.notch.freq_hz;
            bw = cfg.preproc.notch.bw_hz;
            logmsg(logf, '[NOTCH] f0 = %g bw = %g', f0, bw);

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
        [EEG, tags, didLoad] = maybe_load_stage(ST.RESAMPLE, P, subjid, tags, nextTag, logf, EEG);
        if ~didLoad
            targetFs = cfg.preproc.resample.target_hz;
            if EEG.srate ~= targetFs
                logmsg(logf, '[RESAMPLE] %g -> %g Hz', EEG.srate, targetFs);
                EEG = pop_resample(EEG, targetFs);
                EEG = eeg_checkset(EEG);

                tags{end+1} = nextTag;
                save_stage(ST.RESAMPLE, P, subjid, tags, EEG, logf);
            else
                logmsg(logf, '[RESAMPLE] Already at %g Hz (skipping).', targetFs);
            end
        end
    end

    %% -----------
    % REREFERENCE
    % ------------
    if cfg.preproc.reref.enabled
        nextTag = char(string(cfg.preproc.reref.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.REREF, P, subjid, tags, nextTag, logf, EEG);
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
            save_stage(ST.REREF, P, subjid, tags, EEG, logf);
        end
    end

    %% ---------------------------------------------------
    % INITREJ (suggest bad chans + plots + manual interp)
    % ----------------------------------------------------
    if cfg.preproc.initrej.enabled
        nextTag = char(string(cfg.preproc.initrej.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.INITREJ, P, subjid, tags, nextTag, logf, EEG);
        if ~didLoad
            logmsg(logf, '[INITREJ] Suggesting bad channels + manual spherical interpolation.');

            [badChans, reasons, metrics] = suggest_bad_channels(EEG);

            if ~isfield(EEG, 'etc') || isempty(EEG.etc)
                EEG.etc = struct();
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
            else
                logmsg(logf, '[INITREJ] Interpolating channels (spherical): %s', vec2str(interpChans));
                EEG = pop_interp(EEG, interpChans, 'spherical');
                EEG = eeg_checkset(EEG);

                notInterp = setdiff(badChans, interpChans);
                if ~isempty(notInterp)
                    logmsg(logf, '[INITREJ] Suggested but NOT interpolated: %s', vec2str(notInterp));
                    for k = 1:numel(notInterp)
                        ch = notInterp(k);
                        idx = find(badChans == ch, 1);
                        lbl = safe_chan_label(EEG, ch);
                        if ~isempty(idx)
                            logmsg(logf, '  - Ch %d (%s): %s', ch, lbl, reasons{idx});
                        end
                    end
                end
            end

            tags{end+1} = nextTag;
            save_stage(ST.INITREJ, P, subjid, tags, EEG, logf);
        end
    end

    %% ---------------
    % ICA (+ ICLabel)
    % ----------------
    if cfg.preproc.ica.enabled
        nextTag = char(string(cfg.preproc.ica.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.ICA, P, subjid, tags, nextTag, logf, EEG);
        if ~didLoad
            logmsg(logf, '[ICA] method = %s', char(string(cfg.preproc.ica.method)));

            [EEGtrain, segInfo] = make_ica_training_copy(EEG, cfg, logf);

            icatype = char(string(cfg.preproc.ica.method));

            doExtended = false;
            if isfield(cfg.preproc.ica, 'runica') && isfield(cfg.preproc.ica.runica, 'extended')
                doExtended = logical(cfg.preproc.ica.runica.extended);
            end

            if strcmpi(icatype, 'runica') && doExtended
                EEGtrain = pop_runica(EEGtrain, 'icatype', 'runica', 'extended', 1);
            else
                EEGtrain = pop_runica(EEGtrain, 'icatype', icatype);
            end

            EEGtrain = eeg_checkset(EEGtrain);

            logmsg(logf, '[ICA] icatype = %s extended = %d', icatype, doExtended);

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

                    % Save QC packets -> QC folder
                    try
                        save_ic_qc_packets(QC, subjid, EEG, suggestICs);
                        logmsg(logf, '[ICQC] Saved IC QC packets to QC.');
                    catch ME
                        logmsg(logf, '[WARN] IC QC packet generation failed: %s', ME.message);
                    end

                    % Optional tag for ICLabel pass (so filename captures
                    % it)
                    if isfield(cfg.preproc.ica.iclabel, 'tag') && ~isempty(cfg.preproc.ica.iclabel.tag)
                        tags{end+1} = char(string(cfg.preproc.ica.iclabel.tag));
                    end

                    % Numeric mertrics -> LOGS
                    try
                        icMetrics = compute_ic_psd_metrics(EEG, suggestICs);
                        write_ic_metrics_csv(LOGS, subjid, icMetrics);
                        log_ic_metrics(logf, icMetrics);
                        logmsg(logf, '[ICMET] Logged PSD metrics + wrote CSV.');
                    catch ME
                        logmsg(logf, '[WARN] IC PSD metric computation failed: %s', ME.message);

                        logmsg(logf, '[WARN] ICLabel stack trace:\n%s', getReport(ME, 'extended', 'hyperlinks', 'off'));

                        if ~isempty(ME.stack)
                            s0 = ME.stack(1);
                            logmsg(logf, '[WARN] Top stack frame: %s (line %d) | func = %s', s0.file, s0.line, s0.name);
                        end
                        
                    end

                    % Manual reject (default none)
                    removedICs = prompt_ic_reject(suggestICs);

                    if isempty(removedICs)
                        logmsg(logf, '[ICREJ] No ICs removed.');
                    else
                        logmsg(logf, '[ICREJ] Removing ICs: %s', vec2str(removedICs));
                        EEG = pop_subcomp(EEG, removedICs, 0);
                        EEG = eeg_checkset(EEG);
                    end

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

                    logmsg(logf, '[WARN] ICLabel stack trace:\n%s', getReport(ME, 'extended', 'hyperlinks', 'off'));

                    if ~isempty(ME.stack)
                        s0 = ME.stack(1);
                        logmsg(logf, '[WARN] Top stack frame: %s (lind %d) | func = %s', s0.file, s0.line, s0.name);
                    end
                end
            else
                logmsg(logf, '[ICLABEL] disabled.');
            end

            save_stage(ST.ICA, P, subjid, tags, EEG, logf);
        end
    end

    %% -----
    % EPOCH
    % ------
    if cfg.preproc.epoch.enabled
        nextTag = char(string(cfg.preproc.epoch.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.EPOCH, P, subjid, tags, nextTag, logf, EEG);
        if ~didLoad
            tmin = cfg.preproc.epoch.tmin_sec;
            tmax = cfg.preproc.epoch.tmax_sec;

            % ---------------------------------------------
            % Backward compatible event specification
            % ----------------------------------------------
            eventMap = table();

            if isfield(cfg.preproc.epoch, 'events') && ~isempty(cfg.preproc.epoch.events)
                evSpec = cfg.preproc.epoch.events;

                if isstruct(evSpec)
                    nEv = numel(evSpec);
                    code = strings(nEv, 1);
                    condition = strings(nEv, 1);
                    intensity = strings(nEv, 1);

                    for ii = 1:nEv
                        code(ii) = string(evSpec(ii).code);

                        if isfield(evSpec(ii), 'condition') && ~isempty(evSpec(ii).condition)
                            condition(ii) = string(evSpec(ii).condition);
                        else
                            condition(ii) = "unknown";
                        end

                        if isfield(evSpec(ii), 'intensity') && ~isempty(evSpec(ii).intensity)
                            intensity(ii) = string(evSpec(ii).intensity);
                        else
                            intensity(ii) = "unknown";
                        end
                    end

                    eventMap = table(code, condition, intensity);
                else
                    error('cfg.preproc.epoch.events must be a struct array.');
                end

            elseif isfield(cfg.preproc.epoch, 'event_types') && ~isempty(cfg.preproc.epoch.event_types)
                % Legacy mode
                ev = cfg.preproc.epoch.event_types;
                if isstring(ev); ev = cellstr(ev); end
                code = string(ev(:));
                condition = repmat("unknown", size(code), 1);
                intensity = repmat("unknown", size(code), 1);
                eventMap = table(code, condition, intensity);

            else
                error('Epoch config must define either preproc.epoch.events or preproc.epoch.event_types.');
            end

            % Unique event codes used for actual epoching
            ev = cellstr(unique(eventMap.code, 'stable'));

            logmsg(logf, '[EPOCH] events = %s window = [%.3f %.3f] sec', ...
                strjoin(string(ev), ', '), tmin, tmax);

            validate_events_before_epoch(EEG, ev, logf);

            EEG = pop_epoch(EEG, ev, [tmin tmax]);
            EEG = eeg_checkset(EEG);

            % ----------------------------------------------
            % Save event mapping in EEG.etc
            % ----------------------------------------------
            if ~isfield(EEG, 'etc') || isempty(EEG.etc)
                EEG.etc = struct();
            end

            EEG.etc.epoch = struct();
            EEG.etc.event.event_map = eventMap;
            EEG.etc.epoch.tmin_sec = tmin;
            EEG.etc.epoch.tmax_sec = tmax;

            % ----------------------------------------------
            % Per epoch metadata extraction
            % ----------------------------------------------
            nEp = EEG.trials;
            epoch_index = (1:nEp)';
            event_code = strings(nEp, 1);
            condition = strings(nEp, 1);
            intensity = strings(nEp, 1);

            for ep = 1:nEp
                epochEventType = "";

                % EEGLAB stores per-epoch event types in EEG.epoch(ep).eventtype
                if isfield(EEG, 'epoch') && numel(EEG.epoch) >= ep && isfield(EEG.epoch(ep), 'eventtype')
                    et = EEG.epoch(ep).eventtype;

                    if iscell(et)
                        epochEventType = string(et{1});
                    else
                        epochEventType = string(et);
                    end
                end

                event_code(ep) = epochEventType;

                idx = find(eventMap.code == epochEventType, 1, 'first');
                if ~isempty(idx)
                    condition(ep) = "unknown";
                    intensity(ep) = "unknown";
                else
                    condition(ep) = eventMap.condition(idx);
                    intensity(ep) = eventMap.intensity(idx);
                end
            end

            epochMeta = table(epoch_index, event_code, condition, intensity);
            EEG.etc.epoch.epoch_metadata = epochMeta;

            % ----------------------------------------------
            % Save sidecar CSV for downstream joins
            % ----------------------------------------------
            try
                epochCsv = fullfile(ST.EPOCH, sprintf('sub-%03d_epoch_metadata.csv', subjid));
                writetable(epochMeta, epochCsv);
                logmsg(logf, '[EPOCH] Wrote epoch metadata CSV: %s', epochCsv);
            catch ME
                logmsg(logf, '[WARN] Failed to write epoch metadata CSV: %s', ME.message);
            end

            tags{end+1} = nextTag;
            save_stage(ST.EPOCH, P, subjid, tags, EEG, logf);
        end
    end
    %% --------
    % BASELINE
    % ---------
    if cfg.preproc.baseline.enabled
        nextTag = char(string(cfg.preproc.baseline.tag));
        [EEG, tags, didLoad] = maybe_load_stage(ST.BASE, P, subjid, tags, nextTag, logf, EEG);
        if ~didLoad
            win = cfg.preproc.baseline.window_sec;
            logmsg(logf, '[BASELINE] window = [%.3f %.3f] sec', win(1), win(2));

            EEG = pop_rmbase(EEG, win * 1000); % ms
            EEG = eeg_checkset(EEG);

            tags{end+1} = nextTag;
            save_stage(ST.BASE, P, subjid, tags, EEG, logf);
        end
    end

    %% -------------------------------------------------
    % Stage 09: SLOW-ALPHA HILBERT (instantaneous phase)
    % --------------------------------------------------
    % Gate: cfg.preproc.hilbert.enabled (default false; set in JSON)
    %
    % Why here and not in a standalone script:
    %   - cfg.exp.subjects is already resolved by preproc_default, so
    %     subjects_override works automatically with no extra plumbing
    %   - EEG is still in memory from stage 08_base, avoiding a reload
    %   - Resume logic mirrors all other stages: if the .mat already
    %     exists this block is skipped entirely on reruns
    %
    % slow_hz resolution order:
    %   1. cfg.spectral.alpha.slow_hz (preferred; keeps spectral config
    %      as the single source of truth for alpha bands)
    %   2. cfg.preproc.hilbert.slow_hz (fallback; useful when spectral
    %      block is absent from the JSON, default [8 10])
    if isfield(cfg.preproc, 'hilbert') && isfield(cfg.preproc.hilbert, 'enabled') && ...
        logical(cfg.preproc.hilbert.enabled)

        hilbertDir = fullfile(subRoot, '09_hilbert');
        hilbertLogs = fullfile(hilbertDir, 'LOGS');
        ensure_dir(hilbertDir);
        ensure_dir(hilbertLogs);

        matPath = fullfile(hilbertDir, sprintf('sub-%03d_hilbert_phase.mat', subjid));

        if exist(matPath, 'file')
            logmsg(logf, '[HILBERT] Output already exists. Skipping (resume). %s', matPath);
        else
            % Verify EEG is epoched; stage 09 requires trial dimension
            if EEG.trials <= 1
                logmsg(logf, '[HILBERT][WARN] EEG is not epoched (trials = %d).', ...
                    'Stage 09 requires epochs. Ensure epoch + baseline stages ran before Hilbert.', ...
                    EEG.trials);
            else
                % Resolve slow_hz
                slowHz = [8 10]; % ultimate fallback
                if isfield(cfg, 'spectral') && isfield(cfg.spectral, 'alpha') && ...
                    isfield(cfg.spectral.alpha, 'slow_hz')
                    slowHz = cfg.spectral.alpha.slow_hz;
                elseif isfield(cfg.preproc.hilbert, 'slow_hz') && ... 
                    ~isempty(cfg.preproc.hilbert.slow_hz)
                    slowHz = cfg.preproc.hilbert.slow_hz;
                end

                fs = EEG.srate;
                nChan = EEG.nbchan;
                nTr = EEG.trials;
                nTime = EEG.pnts;

                % Stimulus-onset sample index (t closest to 0 ms)
                timesSec = double(EEG.times(:))' / 1000;
                [~, t0idx] = min(abs(timesSec));

                logmsg(logf, '[HILBERT] Band = [%.0f %.0f] Hz nChan = %d nTrials = %d stimOnsetIdx = %d', ...
                    slowHz(1), slowHz(2), nChan, nTr, t0idx);

                % Zero-phase FIR bandpass – operates on a temporary copy
                % EEG is never modified.
                try
                   EEGfilt = pop_eegfiltnew(EEG, slowHz(1), slowHz(2));
                   EEGfilt = eeg_checkset(EEGfilt);
                catch ME
                    logmsg(logf, '[HILBERT][WARN] FIR bandpass failed: %s (skipping stage 09)', ME.message);
                    EEGfilt = []; 
                end

                if ~isempty(EEGfilt)
                    % Hilbert per channel x trial - store phase as single
                    % to keep file size manageable (~4x smaller than double)
                    phaseSlow = zeros(nChan, nTime, nTr, 'single');

                    for t = 1:nTr
                        X = double(EEGfilt.data(:, :, t)); % [nChan x nTime]
                        for ch = 1:nChan
                            z = hilbert(X(ch, :)'); % analytic signal [nTime x 1]
                            phaseSlow(ch, :, t) = single(angle(z));
                        end
                    end

                    % Save - use -v7.3 for arrays > 2 GB
                    phase_slow = phaseSlow;
                    t_ms = EEG.times;
                    stimOnsetIdx = t0idx;
                    chan_labels = {EEG.chanlocs.labels};
                    slow_hz = slowHz;
                    subjid_save = subjid;

                    save(matPath, 'phase_slow', 't_ms', 'stimOnsetIdx', ...
                        'chan_labels', 'slow_hz', 'subjid_save', 'fs', '-v7.3');

                    logmsg(logf, '[HILBERT] Saved -> %s', matPath);
                end
            end
        end
    end

    logmsg(logf, '===== PREPROC DONE sub-%03d =====', subjid);
end

end