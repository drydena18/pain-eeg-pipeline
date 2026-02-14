function write_channel_psd_csv(LOGS, subjid, EEG, chanPSD)
ensure_dir(LOGS);
csvPath = fullfile(LOGS, sprintf('sub-%03d_chan_psd_metrics.csv', subjid));
fid = fopen(csvPath, 'w');
if fid < 0
    warning('Could not write channel PSD CSV: %s', csvPath);
    return;
end

fprintf(fid, 'chan_idx,label,line_ratio,hf_ratio,drift_ratio,alpha_ratio\n');
for ch = 1:EEG.nbchan
    fprintf(fid, '%d,%s,%.6f,%.6f,%.6f,%.6f\n', ...
        ch, safe_chan_label(EEG, ch), ...
        chanPSD.line_ratio(ch), chanPSD.hf_ratio(ch), chanPSD.drift_ratio(ch), chanPSD.alpha_ratio(ch));
end
fclose(fid);
end