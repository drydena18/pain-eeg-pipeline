function out = st_tags(EEG, tags)
%ST_TAGS st_tags(EEG) returns the cumulative stage tags carried in
%   EEG.etc.pp_tags; st_tags(EEG, tags) returns EEG with them set.
if nargin < 2
    out = {};
    if isfield(EEG, 'etc') && isstruct(EEG.etc) && isfield(EEG.etc, 'pp_tags') && ~isempty(EEG.etc.pp_tags)
        out = cellstr(EEG.etc.pp_tags(:)');
    end
else
    if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc), EEG.etc = struct(); end
    EEG.etc.pp_tags = tags;
    out = EEG;
end
end