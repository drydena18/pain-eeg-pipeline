function subs = normalize_subject_ids(subs_raw)
if isnumeric(subs_raw)
    subs = subs_raw(:);
    subs = unique(round(subs));
    subs(subs <= 0) = [];
    return;
end

subs_str = string(subs_raw);
subs = nan(numel(subs_str), 1);

for i = 1:numel(subs_str)
    tok = regexp(char(subs_str(i)), '(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        subs(i) = str2double(tok{1});
    end
end

if any(isnan(subs))
    badIdx = find(isnan(subs), 1, 'first');
    error('preproc_default:BadSubjectId', ...
        'Could not parse subject id from value "%s".', subs_str(badIdx));
end

subs = unique(round(subs(:)));
subs(subs <= 0) = [];
end