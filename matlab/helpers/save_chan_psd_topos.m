function save_chan_psd_topos(LOGS, subjid, EEG, chanPSD)
if ~has_chanlocs(EEG); return; end
save_topo_metric(LOGS, subjid, EEG, chanPSD.line_ratio, 'LINE_RATIO');
save_topo_metric(LOGS, subjid, EEG, chanPSD.hf_ratio, 'HF_RATIO');
save_topo_metric(LOGS, subjid, EEG, chanPSD.drift_ratio, 'DRIFT_RATIO');
save_topo_metric(LOGS, subjid, EEG, chanPSD.alpha_ratio, 'ALPHA_RATIO');
end