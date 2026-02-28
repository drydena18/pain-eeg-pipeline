function fid = spec_open_log(outLogDir, subjid, stem)
spec_ensure_dir(outLogDir);
logPath = fullfile(outLogDir, sprintf('sub-%03d_%s.log', subjid, stem));
fid = fopen(logPath, 'w');
if fid < 0 
    warning('spec_open_log:Fail', 'Could not open log: %s.', logPath);
    fid = 1;
end
end