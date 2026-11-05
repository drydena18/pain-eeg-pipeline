function [rt, cleanup] = st_runtime(cfg, subj, mode)
%ST_RUNTIME Per-participant runtime for the CNED stage hooks: subject
%   number, preproc_core's per-subject folder scheme, LOGS/QC and a log file
%   (LOGS/sub-XXX_pp_<mode>.log, appended).
P  = cfg.P;
id = str2double(regexp(subj, '\d+', 'match', 'once'));
subRoot = fullfile(char(P.RUN_ROOT), sprintf('sub-%03d', id));

ST = struct();
ST.CONCAT   = fullfile(subRoot, '00_concat');
ST.FILTER   = fullfile(subRoot, '01_filter');
ST.NOTCH    = fullfile(subRoot, '02_notch');
ST.RESAMPLE = fullfile(subRoot, '03_resample');
ST.REREF    = fullfile(subRoot, '04_reref');
ST.INITREJ  = fullfile(subRoot, '05_initrej');
ST.ICA      = fullfile(subRoot, '06_ica');
ST.EPOCH    = fullfile(subRoot, '07_epoch');
ST.BASE     = fullfile(subRoot, '08_base');
ST.HILBERT  = fullfile(subRoot, '09_hilbert');
LOGS = fullfile(subRoot, 'LOGS');
QC   = fullfile(subRoot, 'QC');

ensure_dir(subRoot); ensure_dir(LOGS); ensure_dir(QC);
f = fieldnames(ST);
for k = 1:numel(f), ensure_dir(ST.(f{k})); end

logPath = fullfile(LOGS, sprintf('sub-%03d_pp_%s.log', id, mode));
logf = fopen(logPath, 'a');
if logf < 0
    warning('st_runtime:LogOpenFail', 'Could not open log file: %s', logPath);
    logf = 1;
end
cleanup = onCleanup(@() safeClose(logf));

logmsg(logf, '===== PP %s START sub-%03d (%s) =====', upper(mode), id, string(cfg.exp_id));
rt = struct('subjid', id, 'P', P, 'subRoot', subRoot, 'ST', ST, ...
            'LOGS', LOGS, 'QC', QC, 'logf', logf, 'resume', true);
end