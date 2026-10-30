function [EEG, qc] = st_stage(EEG, cfg, stKey, tag, computeFn)
%ST_STAGE preproc_core's stage pattern. If cfg.rt.resume and the stage file
%   exists, load it (maybe_load_stage); otherwise compute, append the tag and
%   save with save_stage - same folder, same file name.
%   computeFn: [EEG, qc] = computeFn(EEG). Setting qc.no_save = true skips
%   the tag and save (resample when already at the target rate).
rt   = cfg.rt;
tags = st_tags(EEG);
tag  = char(string(tag));

didLoad = false;
if rt.resume
    [E2, tags2, didLoad] = maybe_load_stage(rt.ST.(stKey), rt.P, rt.subjid, tags, tag, rt.logf, EEG);
end

if didLoad
    EEG = st_tags(E2, tags2);
    tags = tags2;
    qc = struct('loaded_existing', true);
else
    [EEG, qc] = computeFn(EEG);
    if isempty(qc), qc = struct(); end
    qc.loaded_existing = false;
    if isfield(qc, 'no_save') && qc.no_save
        EEG = st_tags(EEG, tags);
        return;
    end
    tags{end+1} = tag;
    EEG = st_tags(EEG, tags);
    save_stage(rt.ST.(stKey), rt.P, rt.subjid, tags, EEG, rt.logf);
end
qc.tag   = tag;
qc.saved = fullfile(rt.ST.(stKey), char(rt.P.NAMING.fname(rt.subjid, tags, [])));
end