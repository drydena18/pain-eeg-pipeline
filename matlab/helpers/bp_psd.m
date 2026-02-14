function p = bp_psd(f, pxx, band)
m = (f >= band(1) & f < band(2));
if ~any(m), p = 0; return; end
p = trapz(f(m), pxx(m));
end