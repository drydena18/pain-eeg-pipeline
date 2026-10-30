function pp_quarantine(file)
%PP_QUARANTINE Rename instead of delete, so manual work is never lost.
if ~isfile(file), return; end
[d, n, e] = fileparts(file);
dst = fullfile(d, sprintf('%s_stale_%s%s', n, char(datetime('now'), 'yyyyMMdd''T''HHmmss'), e));
movefile(file, dst);
warning('pp:stale', 'Moved %s -> %s (made for a previous ICA run).', [n e], dst);
end