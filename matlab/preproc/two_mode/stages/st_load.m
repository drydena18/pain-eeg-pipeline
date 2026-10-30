function EEG = st_load(cfg, ~)
%ST_LOAD Stage 00 of preproc_core: multi-session concatenation (saved to and
%   resumable from 00_concat) or single-session load, then BIDS channel
%   labels, montage, ELP channel locations and the standard-label guard.
rt = cfg.rt; ecfg = cfg.exp; P = rt.P; subjid = rt.subjid; logf = rt.logf;
tags = {};

logmsg(logf, 'Experiment: %s | Out Prefix: %s | Output Root: %s', ...
    string(ecfg.exp.id), string(ecfg.exp.out_prefix), rt.subRoot);

doConcat = isfield(ecfg.preproc, 'concat') && isfield(ecfg.preproc.concat, 'enabled') && ...
           logical(ecfg.preproc.concat.enabled);

if doConcat
    nSess = ecfg.preproc.concat.n_sessions;
    concatTag = 'concat';
    logmsg(logf, '[CONCAT] Multi-session mode: %d sessions.', nSess);

    didLoad = false;
    if rt.resume
        [EEG, tags, didLoad] = maybe_load_stage(rt.ST.CONCAT, P, subjid, tags, concatTag, logf);
    end

    if ~didLoad
        [sessPaths, sessSidecars] = resolve_raw_session_files(P, ecfg, subjid, nSess, logf);
        sessEEGs = cell(nSess, 1);
        for s = 1:nSess
            logmsg(logf, '[CONCAT] Loading session %d: %s', s, sessPaths{s});
            Etmp = load_raw_eeg(sessPaths{s}, logf);
            Etmp = eeg_checkset(Etmp);
            Etmp = normalize_chan_labels(Etmp);

            sidecar = sessSidecars{s};
            if ~isempty(sidecar.channels_tsv)
                try
                    Etmp = apply_bids_channel_labels(Etmp, sidecar.channels_tsv, logf);
                    Etmp = eeg_checkset(Etmp);
                catch ME
                    logmsg(logf, '[CONCAT][WARN] channels.tsv apply failed sess %d: %s', s, ME.message);
                end
            end

            Etmp = local_montage(Etmp, P, ecfg, logf, subjid, rt.LOGS);
            Etmp = local_chanlocs(Etmp, P, ecfg, logf, sprintf('[CONCAT] sess %d', s), false);

            logmsg(logf, '[CONCAT] Session %d: %d chans  %.1f s  srate=%.0f Hz', ...
                s, Etmp.nbchan, Etmp.xmax, Etmp.srate);
            sessEEGs{s} = Etmp;
        end

        logmsg(logf, '[CONCAT] Merging %d sessions with pop_mergeset...', nSess);
        EEG = sessEEGs{1};
        for s = 2:nSess
            EEG = pop_mergeset(EEG, sessEEGs{s});
        end
        EEG = eeg_checkset(EEG);
        logmsg(logf, '[CONCAT] Merged: %d chans  %.1f s  %d events', EEG.nbchan, EEG.xmax, length(EEG.event));

        EEG.etc.concat = struct();
        EEG.etc.concat.n_sessions       = nSess;
        EEG.etc.concat.session_files    = sessPaths;
        EEG.etc.concat.session_sidecars = sessSidecars;

        tags{end+1} = concatTag;
        EEG = st_tags(EEG, tags);
        save_stage(rt.ST.CONCAT, P, subjid, tags, EEG, logf);
    end
else
    rawPath = resolve_raw_file(P, ecfg, subjid);
    if isempty(rawPath) || ~exist(rawPath, 'file')
        error('st_load:RawMissing', 'Raw file not found for sub-%03d.', subjid);
    end
    rawPath = char(string(rawPath));
    logmsg(logf, 'Raw File: %s', rawPath);

    EEG = load_raw_eeg(rawPath, logf);
    EEG = eeg_checkset(EEG);
    EEG = normalize_chan_labels(EEG);
    EEG = local_montage(EEG, P, ecfg, logf, subjid, rt.LOGS);
    EEG = local_chanlocs(EEG, P, ecfg, logf, '', true);
end

% Guard: no raw BioSemi / auxiliary labels past this point (errors stop this participant)
check_standard_chan_labels(EEG, logf, sprintf('sub-%03d:', subjid));
EEG = st_tags(EEG, tags);
end

% -------------------------------------------------------------------------
function EEG = local_montage(EEG, P, ecfg, logf, subjid, LOGS)
if isfield(ecfg.exp, 'montage') && isfield(ecfg.exp.montage, 'enabled') && ecfg.exp.montage.enabled
    EEG = apply_montage_biosemi_from_csv(P, ecfg, EEG, logf, subjid, LOGS);
    EEG = eeg_checkset(EEG);
end
end

function EEG = local_chanlocs(EEG, P, ecfg, logf, ctx, warnMissing)
if ~(isfield(ecfg.exp, 'channel_locs') && isfield(ecfg.exp.channel_locs, 'use_elp') && ...
        ecfg.exp.channel_locs.use_elp)
    return;
end
elpPath = "";
if isfield(P, 'CORE') && isfield(P.CORE, 'ELP_FILE'), elpPath = string(P.CORE.ELP_FILE); end
if strlength(elpPath) > 0 && exist(elpPath, 'file')
    try
        EEG = pop_chanedit(EEG, 'lookup', char(elpPath));
        EEG = eeg_checkset(EEG);
        logmsg(logf, '%s Loaded chanlocs from ELP: %s', ctx, elpPath);
    catch ME
        logmsg(logf, '[WARN]%s Failed to load chanlocs: %s', ctx, ME.message);
    end
elseif warnMissing
    logmsg(logf, '[WARN] ELP file missing: %s', elpPath);
end
end