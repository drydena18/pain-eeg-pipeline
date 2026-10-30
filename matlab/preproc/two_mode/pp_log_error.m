function pp_log_error(cfg, subj, mode, ME)
%PP_LOG_ERROR Append the full stack trace to logs/pp_errors.log.
P = pp_paths(cfg, subj);
fid = fopen(P.log, 'a');
if fid < 0, return; end
fprintf(fid, '[%s] %s  Mode %s\n%s\n\n', char(datetime('now'), 'yyyy-MM-dd HH:mm:ss'), ...
    subj, mode, ME.getReport('extended', 'hyperlinks', 'off'));
fclose(fid);
end