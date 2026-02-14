function logmsg(fid)
if fid ~= 1 && fid > 0
    fclose(fid);
end
end