function [EEG, qc] = st_notch(EEG, cfg, ~)
%ST_NOTCH Stage 02 (preproc_core NOTCH block, bandstop via revfilt) -> 02_notch.
n = cfg.exp.preproc.notch;
if ~logical(n.enabled), qc = struct('disabled', true); return; end
[EEG, qc] = st_stage(EEG, cfg, 'NOTCH', n.tag, @compute);

    function [E, q] = compute(E)
        f0 = n.freq_hz; bw = n.bw_hz;
        logmsg(cfg.rt.logf, '[NOTCH] f0 = %g bw = %g', f0, bw);
        E = pop_eegfiltnew(E, f0 - bw, f0 + bw, [], 1);   % revfilt = 1 -> bandstop
        E = eeg_checkset(E);
        q = struct('freq_hz', f0, 'bw_hz', bw);
    end
end