function s = vec2str(v)
if isempty(v), s = '[]'; return, end
v = v(:)';
s = ['[' sprintf('%d ', v) ']'];
s = strrep(s, ' ]', ']');
end