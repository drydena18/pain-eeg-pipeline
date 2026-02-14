function mustHave(s, field, msg)
if ~isfield(s, field)
    error('preproc_default:MissingField', '%s', msg);
end
end