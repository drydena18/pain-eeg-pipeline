function make_initrej_plots(LOGS, subjid, EEG, metrics, badChans)
ensure_dir(LOGS);

save_hist(LOGS, subjid, metrics.chan_std, 'chan_std');
save_hist(LOGS, subjid, metrics.chan_rms, 'chan_rms');

save_bar(LOGS, subjid, metrics.chan_std, badChans, 'chan_std');
save_bar(LOGS, subjid, metrics.chan_rms, badChans, 'chan_rms');

if has_chanlocs(EEG)
    save_topo_metric(LOGS, subjid, EEG, metrics.chan_std, 'STD');
    save_topo_metric(LOGS, subjid, EEG, metrics.chan_rms, 'RMS');
end

save_channel_psd_overview(LOGS, subjid, EEG);

if ~isempty(badChans)
    save_channel_psd_badchans(LOGS, subjid, EEG, badChans);
end

% Label index map
fid = fopen(fullfile(LOGS, sprintf('sub-%03d_chalabels.txt', subjid)), 'w');
if fid > 0
    for k = 1:EEG.nbchan
        fprintf(fid, '%d\t%s\n', k, safe_chan_label(EEG, k));
    end
    fclose(fid);
end
end