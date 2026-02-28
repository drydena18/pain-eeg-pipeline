function labels = spec_get_chanlabels(EEG)
labels = strings(EEG.nbchan, 1);
if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs, 'labels')
    for i = 1:EEG.nbchan
        labels(i) = string(EEG.chanlocs(i).labels);
    end
else
    for i = 1:EEG.nbchan
        labels(i) = "Ch" + string(i);
    end
end
end