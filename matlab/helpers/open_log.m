function [fid, logPath, cleanupObj] = open_log(LOGS, subjid, stem)
ensure_dir(LOGS);
if nargin < 3 || isempty(stem), stem = 'log'; end
logPath = fullfile(LOGS, sprintf('sub-%03d_%s.log', subjid, stem));
fid = fopen(logPath, 'w');
if fid < 0
    warning('open_log:Fail', 'Could not open log file: %s (using stdout)', logPath);
    fid = 1;
end
cleanupObj = onCleanup(@() safeClose(fid));
end