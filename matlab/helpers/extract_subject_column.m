function subs_raw = extract_subject_column(T)
vn = lower(string(T.Properties.VariableNames));

if any(vn == 'subjid')
    subs_raw = T.(T.Properties.VariableNames{find(vn == 'subjid', 1)});
elseif any(vn == 'participant_id')
    subs_raw = T.(T.Properties.VariableNames{find(vn == 'participant_id', 1)});
elseif any(vn == 'subject')
    subs_raw = T.(T.Properties.VariableNames{find(vn == 'subject', 1)});
else
    error('preproc_default:BadParticipantsTSV', ...
        'participants.tsv needs a column like subjid / participant_id / subject.');
end
end
