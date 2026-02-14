function EEG = ensure_etc_path(EEG)
if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
end