function [EEG, qc] = st_resample(EEG, cfg, ~)
%ST_RESAMPLE Stage 03 (preproc_core RESAMPLE block) -> 03_resample.
%   Like preproc_core, nothing is saved (and no tag added) if the data are
%   already at the target rate.
r = cfg.exp.preproc.resample;
if ~logical(r.enabled), qc = struct('disabled', true); return; end
[EEG, qc] = st_stage(EEG, cfg, 'RESAMPLE', r.tag, @compute);

    function [E, q] = compute(E)
        targetFs = r.target_hz;
        q = struct('srate_in', E.srate, 'target_hz', targetFs);
        if E.srate ~= targetFs
            logmsg(cfg.rt.logf, '[RESAMPLE] %g -> %g Hz', E.srate, targetFs);
            E = pop_resample(E, targetFs);
            E = eeg_checkset(E);
        else
            logmsg(cfg.rt.logf, '[RESAMPLE] Already at %g Hz (skipping).', targetFs);
            q.no_save = true;
        end
    end
end