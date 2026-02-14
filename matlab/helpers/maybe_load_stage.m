function [EEG, tags, didLoad] = maybe_load_stage(stageDir, P, subjid, tags, nextTag, logf)
didLoad = false;
tags2 = tags;
if ~isempty(nextTag)
    tags2{end+1} = nextTag;
end

fname = P.NAMING.fname(subjid, tags2, []);
fpath = fullfile(stageDir, fname);

if exist(fpath, 'file')
    logmsg(logf, '[SKIP] Loading existing stage file: %s', fpath);
    EEG = pop_loadset('filename', fname, 'filepath', stageDir);
    EEG = eeg_checkset;
    tags = tags2;
    didLoad = true;
end
end