function [EEG, qc] = st_baseline(EEG, cfg, ~)
%ST_BASELINE Stage 08 (preproc_core BASELINE block) -> 08_base. This is the
%   file spectral (input_stage "08_base") and source localization read.
b = cfg.exp.preproc.baseline;
if ~logical(b.enabled), qc = struct('disabled', true); return; end
[EEG, qc] = st_stage(EEG, cfg, 'BASE', b.tag, @compute);

    function [E, q] = compute(E)
        win = b.window_sec;
        logmsg(cfg.rt.logf, '[BASELINE] window = [%.3f %.3f] sec', win(1), win(2));
        E = pop_rmbase(E, win * 1000);   % ms
        E = eeg_checkset(E);
        q = struct('window_sec', win(:)');
    end
end