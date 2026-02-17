function logmsg(fid, fmt, varargin)
ts = datetime("now", "Format", 'yyyy-MM-dd HH:mm:ss');
msg = sprintf(fmt, varargin{:});
fprintf(fid, '[%s] %s\n', char(ts), msg);
if fid ~= 1
    fprintf(1, '[%s] %s\n', char(ts), msg);
end
end