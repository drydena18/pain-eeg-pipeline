function [EEG, suggestICs, qc] = st_icqc(EEG, cfg, ~)
%ST_ICQC Mode 2a prep: the automated half of preproc_core's ICLabel block.
%   iclabel -> iclabel_suggest_reject (JSON thresholds) -> QC packets
%   (QC/sub-XXX_icqc) -> IC PSD metrics CSV (LOGS) -> first-N topo grid.
rt = cfg.rt; logf = rt.logf; subjid = rt.subjid;
ica = cfg.exp.preproc.ica;
if ~(isfield(ica, 'iclabel') && isfield(ica.iclabel, 'enabled') && ica.iclabel.enabled)
    error('st_icqc:ICLabelOff', 'Mode 2a needs preproc.ica.iclabel.enabled = true.');
end

EEG = iclabel(EEG);
EEG = eeg_checkset(EEG);

thr = struct();
if isfield(ica.iclabel, 'thresholds'), thr = ica.iclabel.thresholds; end
[suggestICs, icReasons] = iclabel_suggest_reject(EEG, thr);
logmsg(logf, '[ICLABEL] Suggested ICs: %s', vec2str(suggestICs));
for k = 1:numel(icReasons), logmsg(logf, '  - %s', icReasons{k}); end

try
    save_ic_qc_packets(rt.QC, subjid, EEG, suggestICs);
    logmsg(logf, '[ICQC] Saved IC QC packets to QC.');
catch ME
    logmsg(logf, '[WARN] IC QC packet generation failed: %s', ME.message);
end
try
    icMetrics = compute_ic_psd_metrics(EEG, suggestICs);
    write_ic_metrics_csv(rt.LOGS, subjid, icMetrics);
    log_ic_metrics(logf, icMetrics);
    logmsg(logf, '[ICMET] Logged PSD metrics + wrote CSV.');
catch ME
    logmsg(logf, '[WARN] IC PSD metric computation failed: %s', ME.message);
end
if isfield(ica, 'grid') && isfield(ica.grid, 'enabled') && ica.grid.enabled
    try
        save_ic_topo_grid(rt.QC, subjid, EEG, ica.grid.n_ics);
        logmsg(logf, '[ICQC] Saved first-%d IC topomap grid to QC.', ica.grid.n_ics);
    catch ME
        logmsg(logf, '[WARN] IC topomap grid generation failed: %s', ME.message);
    end
end

qc = struct();
qc.reasons = icReasons;
end