function [EEG, qc] = st_filter(EEG, cfg, ~)
%ST_FILTER Stage 01 (preproc_core FILTER block) -> 01_filter.
f = cfg.exp.preproc.filter;
if ~logical(f.enabled), qc = struct('disabled', true); return; end
[EEG, qc] = st_stage(EEG, cfg, 'FILTER', f.tag, @compute);

    function [E, q] = compute(E)
        logmsg(cfg.rt.logf, '[FILTER] %s hp = %.3f lp = %.3f', string(f.type), f.highpass_hz, f.lowpass_hz);
        E = pop_eegfiltnew(E, f.highpass_hz, f.lowpass_hz);
        E = eeg_checkset(E);
        q = struct('highpass_hz', f.highpass_hz, 'lowpass_hz', f.lowpass_hz);
    end
end