function pp_write_json(file, s)
%PP_WRITE_JSON Write to a temp file, then rename. A crash mid-write leaves
%   only the .tmp file, never a half-written marker.
txt = jsonencode(pp_sanitize(s), 'PrettyPrint', true);
tmp = [file '.tmp'];
fid = fopen(tmp, 'w');
if fid < 0, error('pp:io', 'Cannot write %s', tmp); end
fwrite(fid, unicode2native(txt, 'UTF-8'), 'uint8');
fclose(fid);
movefile(tmp, file, 'f');
end