function [EEG, qc] = st_epoch(EEG, cfg, ~)
%ST_EPOCH Stage 07 (preproc_core EPOCH block, V2.2.0 event mapping) ->
%   07_epoch, plus sub-XXX_epoch_metadata.csv.
e = cfg.exp.preproc.epoch;
if ~logical(e.enabled), qc = struct('disabled', true); return; end
[EEG, qc] = st_stage(EEG, cfg, 'EPOCH', e.tag, @compute);

    function [E, q] = compute(E)
        rt = cfg.rt; logf = rt.logf; subjid = rt.subjid;
        tmin = e.tmin_sec;
        tmax = e.tmax_sec;

        % ---- Event specification (struct array, or legacy event_types) ----
        if isfield(e, 'events') && ~isempty(e.events)
            evSpec = e.events;
            if ~isstruct(evSpec), error('st_epoch:Events', 'preproc.epoch.events must be a struct array.'); end
            nEv = numel(evSpec);
            code = strings(nEv, 1); condition = strings(nEv, 1); intensity = strings(nEv, 1);
            for ii = 1:nEv
                code(ii) = string(evSpec(ii).code);
                condition(ii) = "unknown"; intensity(ii) = "unknown";
                if isfield(evSpec(ii), 'condition') && ~isempty(evSpec(ii).condition), condition(ii) = string(evSpec(ii).condition); end
                if isfield(evSpec(ii), 'intensity') && ~isempty(evSpec(ii).intensity), intensity(ii) = string(evSpec(ii).intensity); end
            end
            eventMap = table(code, condition, intensity);
        elseif isfield(e, 'event_types') && ~isempty(e.event_types)
            evl = e.event_types;
            if isstring(evl), evl = cellstr(evl); end
            code = string(evl(:));
            condition = repmat("unknown", size(code), 1);
            intensity = repmat("unknown", size(code), 1);
            eventMap = table(code, condition, intensity);
        else
            error('st_epoch:Events', 'Epoch config must define preproc.epoch.events or preproc.epoch.event_types.');
        end

        ev = cellstr(unique(eventMap.code, 'stable'));
        logmsg(logf, '[EPOCH] events = %s window = [%.3f %.3f] sec', strjoin(string(ev), ', '), tmin, tmax);
        validate_events_before_epoch(E, ev, logf);

        E = pop_epoch(E, ev, [tmin tmax]);
        E = eeg_checkset(E);

        if ~isfield(E, 'etc') || isempty(E.etc), E.etc = struct(); end
        E.etc.epoch = struct();
        E.etc.epoch.event_map = eventMap;
        E.etc.epoch.tmin_sec  = tmin;
        E.etc.epoch.tmax_sec  = tmax;

        % ---- Per-epoch metadata from the TIME-LOCKING event (latency 0) ----
        nEp = E.trials;
        epoch_index = (1:nEp)';
        event_code = strings(nEp, 1);
        econd = strings(nEp, 1);
        eint  = strings(nEp, 1);
        for ep = 1:nEp
            epochEventType = "";
            if isfield(E, 'epoch') && numel(E.epoch) >= ep && isfield(E.epoch(ep), 'eventtype')
                et = E.epoch(ep).eventtype;
                if ~iscell(et), et = {et}; end
                etStr = strings(1, numel(et));
                for qq = 1:numel(et), etStr(qq) = string(et{qq}); end
                pick = 1;
                if isfield(E.epoch(ep), 'eventlatency') && ~isempty(E.epoch(ep).eventlatency)
                    lat = E.epoch(ep).eventlatency;
                    if ~iscell(lat), lat = num2cell(lat); end
                    latv = cellfun(@(v) double(v(1)), lat);
                    cand = find(abs(latv) < 1e-6 & ismember(etStr, eventMap.code), 1);
                    if isempty(cand), [~, cand] = min(abs(latv)); end
                    pick = cand;
                end
                epochEventType = etStr(pick);
            end
            event_code(ep) = epochEventType;
            idx = find(eventMap.code == epochEventType, 1, 'first');
            if isempty(idx)
                econd(ep) = "unknown"; eint(ep) = "unknown";
            else
                econd(ep) = eventMap.condition(idx); eint(ep) = eventMap.intensity(idx);
            end
        end
        condition = econd; intensity = eint;
        epochMeta = table(epoch_index, event_code, condition, intensity);
        E.etc.epoch.epoch_metadata = epochMeta;

        nUnknown = sum(condition == "unknown");
        if nUnknown > 0
            logmsg(logf, '[EPOCH][WARN] %d / %d epochs could not be mapped to a configured event code.', nUnknown, nEp);
        end

        try
            epochCsv = fullfile(rt.ST.EPOCH, sprintf('sub-%03d_epoch_metadata.csv', subjid));
            writetable(epochMeta, epochCsv);
            logmsg(logf, '[EPOCH] Wrote epoch metadata CSV: %s', epochCsv);
        catch ME
            logmsg(logf, '[WARN] Failed to write epoch metadata CSV: %s', ME.message);
        end

        q = struct('n_epochs', nEp, 'n_unknown', nUnknown, 'tmin_sec', tmin, 'tmax_sec', tmax);
    end
end