function [EEG, qc] = st_reref(EEG, cfg, ~)
%ST_REREF Stage 04 (preproc_core REREFERENCE block) -> 04_reref.
r = cfg.exp.preproc.reref;
if ~logical(r.enabled), qc = struct('disabled', true); return; end
[EEG, qc] = st_stage(EEG, cfg, 'REREF', r.tag, @compute);

    function [E, q] = compute(E)
        mode = lower(string(r.mode));
        logmsg(cfg.rt.logf, '[REREF] mode = %s', mode);
        if mode == "average"
            E = pop_reref(E, []);
        elseif mode == "channels"
            if isempty(r.channels)
                error('st_reref:NoChannels', 'reref.mode = "channels" but preproc.reref.channels is empty.');
            end
            E = pop_reref(E, r.channels);
        else
            error('st_reref:BadMode', 'Unsupported reref.mode: %s', mode);
        end
        E = eeg_checkset(E);
        q = struct('mode', char(mode));
    end
end