function spec_safe_close(fid)
if fid ~= 1 && fid > 0
    fclose(fid);
end
end