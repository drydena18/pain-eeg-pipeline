function tf = has_chanlocs(EEG)
tf = isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs(1), 'X');
end