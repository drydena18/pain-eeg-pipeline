function lbl = safe_chan_label(EEG, ch)
lbl = '';
try
    if isfield(EEG, 'chanlocs') && numel(EEG.chanlocs) >= ch && isfield(EEG.chanlocs(ch), 'labels')
        lbl = EEG.chanlocs(ch).labels;
    end
catch
end
if isempty(lbl)
    lbl = sprintf('Ch%d', ch);
end
end