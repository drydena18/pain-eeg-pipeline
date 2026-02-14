function EEG = normalize_chan_labels(EEG)
if ~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs), return; end
for i = 1:numel(EEG.chanlocs)
    if isfield(EEG.chanlocs(i), 'labels') && ~isempty(EEG.chanlocs(i).labels)
        lbl = string(EEG.chanlocs(i).labels);
        lbl = strtrim(lbl);
        lbl = replace(lbl, " ", "");
        lbl = replace(lbl, "-", "");
        EEG.chanlocs(i).labels = char(lbl);
    end
end
end