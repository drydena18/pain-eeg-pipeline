function save_stage(stageDir, P, subjid, tags, EEG, logf)

stageDir = char(string(stageDir));
ensure_dir(stageDir);

fname = char(string(P.NAMING.fname(subjid, tags, [])));
outPath = fullfile(stageDir, fname);

logmsg(logf, '[SAVE] %s', outPath);

pop_saveset(EEG, 'filename', fname, 'filepath', stageDir);

end