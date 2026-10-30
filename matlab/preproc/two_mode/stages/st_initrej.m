function [EEG, qc] = st_initrej(EEG, cfg, ~)
%ST_INITREJ Stage 05 (preproc_core INITREJ block) -> 05_initrej, unattended.
%   Same bad-channel suggestions, channel PSD metrics/CSV/topos and QC plots.
%   The channel scroll and prompt_channel_interp are replaced by
%   cfg.mode1.interp_policy:
%     'suggested'  interpolate every channel suggest_bad_channels flags
%     'none'       interpolate nothing (the prompt's Enter default)
%   The decision is stored in EEG.etc.initrej (interp_chans, interp_labels,
%   decision) so the ICA rank can be computed later from the stage file.
ir = cfg.exp.preproc.initrej;
if ~logical(ir.enabled)
    qc = struct('disabled', true);
    qc.bad_labels = {}; qc.interp_labels = {};
    return;
end
policy = lower(char(cfg.mode1.interp_policy));
[EEG, qc] = st_stage(EEG, cfg, 'INITREJ', ir.tag, @compute);

if qc.loaded_existing
    % Reused 05_initrej file (e.g. an earlier manual preproc_core run).
    ini = struct();
    if isfield(EEG.etc, 'initrej'), ini = EEG.etc.initrej; end
    qc.bad_labels = {};
    if isfield(ini, 'suggested_badchans')
        qc.bad_labels = local_labels(EEG, ini.suggested_badchans);
    end
    if isfield(ini, 'interp_labels')
        qc.interp_labels = cellstr(ini.interp_labels);
        if isfield(ini, 'decision'), qc.decision = ini.decision; end
    else
        qc.interp_labels  = {};
        qc.interp_unknown = true;    % preproc_core never stored the manual choice
        qc.decision       = 'manual (preproc_core, not recorded)';
    end
end

    function [E, q] = compute(E)
        rt = cfg.rt; logf = rt.logf; subjid = rt.subjid;
        logmsg(logf, '[INITREJ] Suggesting bad channels; unattended interpolation policy = %s.', policy);

        [badChans, reasons, metrics] = suggest_bad_channels(E);
        if ~isfield(E, 'etc') || isempty(E.etc), E.etc = struct(); end
        E.etc.initrej = struct();
        E.etc.initrej.suggested_badchans = badChans;
        E.etc.initrej.reasons = reasons;
        E.etc.initrej.metrics = metrics;

        try
            chanPSD = compute_channel_psd_metrics(E);
            E.etc.initrej.chan_psd = chanPSD;
            write_channel_psd_csv(rt.LOGS, subjid, E, chanPSD);
            save_chan_psd_topos(rt.LOGS, subjid, E, chanPSD);
            logmsg(logf, '[INITREJ] Channel PSD metrics saved (CSV + topos).');
        catch ME
            logmsg(logf, '[WARN] Channel PSD metrics failed: %s', ME.message);
        end
        try
            make_initrej_plots(rt.LOGS, subjid, E, metrics, badChans);
            logmsg(logf, '[INITREJ] Saved QC plots to LOGS.');
        catch ME
            logmsg(logf, '[WARN] INITREJ plotting failed: %s', ME.message);
        end

        if isempty(badChans)
            logmsg(logf, '[INITREJ] No channel suggested.');
        else
            logmsg(logf, '[INITREJ] Suggested bad channels: %s', vec2str(badChans));
            for k = 1:numel(badChans)
                logmsg(logf, '  - Ch %d (%s): %s', badChans(k), safe_chan_label(E, badChans(k)), reasons{k});
            end
        end

        switch policy
            case 'suggested', interpChans = unique(badChans(:))';
            case 'none',      interpChans = [];
            otherwise, error('st_initrej:Policy', 'Unknown interp_policy "%s".', policy);
        end
        interpLabels = local_labels(E, interpChans);
        E.etc.initrej.interp_chans  = interpChans;
        E.etc.initrej.interp_labels = interpLabels;
        E.etc.initrej.decision      = ['auto:' policy];

        if isempty(interpChans)
            logmsg(logf, '[INITREJ] No channels interpolated (policy = %s).', policy);
        else
            logmsg(logf, '[INITREJ] Interpolating channels (spherical): %s', vec2str(interpChans));
            E = pop_interp(E, interpChans, 'spherical');
            E = eeg_checkset(E);
        end

        q = struct('policy', policy);
        q.bad_labels    = local_labels(E, badChans);
        q.interp_labels = interpLabels;
    end
end

function c = local_labels(EEG, idx)
c = arrayfun(@(ch) char(string(safe_chan_label(EEG, ch))), idx(:)', 'UniformOutput', false);
end