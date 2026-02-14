function ensure_dir(d)
d = char(string(d));
if ~exist(d, 'dir')
    mkdir(d);
end
end