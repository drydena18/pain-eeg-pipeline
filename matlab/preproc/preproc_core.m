function preproc_core(P, cfg)
% PREPROC_CORE Execute preprocessing according to cfg (from JSON)
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

    rawPath = char(string(rawPath));
    assert(exist(rawPath, 'file') == 2, 'Raw file missing: %s', rawPath);

    logmsg(logf, 'RawPath class = %s | isstring = %d', class(rawPath), isstring(rawPath));

    % -------------------------------------------------
    % Import raw (.eeg/.bdf) -> EEG struct (no saving)
    % -------------------------------------------------
    EEG = pop_biosig(rawPath);
    EEG = eeg_checkset(EEG);

    EEG = normalize_chan_labels(EEG);

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
    end

    % Deterministic seed per subject (repro ICA runs)
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

    logmsg(logf, '===== PREPROC DONE sub-%03d =====', subjid);
end

end