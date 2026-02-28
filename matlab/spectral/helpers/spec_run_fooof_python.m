function fooofOut = spec_run_fooof_python(f, gaPxx, foo, outTmp, subjid, logf)
% Runs python fooof_bridge.py on GA PSD per trial
% Inputs:
%   f       [1 x nFreq]
%   gaPxx   [nFreq x nTrial]
%
% Writes:
%   tmp/sub-###_fooof_freqs.csv
%   tmp/sub-###_fooof_psd.csv
%   tmp/sub-###_fooof_cfg.json
%   tmp/sub-###_fooof_out.json

spec_ensure_dir(outTmp);

nF = numel(f);
[nF2, nTr] = size(gaPxx);
if nF2 ~= nF
    error('spec_run_fooof_python:Shape', 'gaPxx must be [nFreq x nTrial].');
end

freqPath = fullfile(outTmp, sprintf('sub-%03d_fooof_freqs.csv', subjid));
psdPath  = fullfile(outTmp, sprintf('sub-%03d_fooof_psd.csv', subjid));
cfgPath  = fullfile(outTmp, sprintf('sub-%03d_fooof_cfg.json', subjid));
outPath  = fullfile(outTmp, sprintf('sub-%03d_fooof_out.json', subjid));

% Write inputs
writematrix(f(:), freqPath);
writematrix(gaPxx', psdPath); % trials x freqs

% Write cfg json for python
pcfg = struct();
pcfg.fmin_hz = foo.min_hz;
pcfg.fmax_hz = foo.max_hz;
pcfg.peak_width_limits = foo.peak_width_limits;
pcfg.max_n_peaks = foo.max_n_peaks;
pcfg.min_peak_height = foo.min_peak_height;
pcfg.aperiodic_mode = char(string(foo.aperiodic_mode));
pcfg.alpha_band_hz = foo.alpha_band_hz;
pcfg.verbose = isfield(foo, 'verbose') && logical(foo.verbose);

txt = jsonencode(pcfg);
fid = fopen(cfgPath, 'w'); fwrite(fid, txt); fclose(fid);

% Resolve python exe + script path
pyexe = "python";
if isfield(foo, 'python_exe') && strlength(string(foo.pythin_exe)) > 0
    pyexe = string(foo.python_exe);
end

script = string(foo.script_path);
if strlength(script) == 0
    error('spec_run_fooof_python:MissingScript', 'cfg.spectral.fooof.script_path is required.');
end

cmd = sprintf('"%s" "%s" --freq "%s" --psd "%s" --cfg "%s" --out "%s"', ...
    pyexe, script, freqPath, psdPath, cfgPath, outPath);

spec_logmsg(logf, '[FOOOF] CMD: %s', cmd);
[status, out] = system(cmd);
spec_logmsg(logf, '[FOOOF] status = %d', status);
if ~isempty(out), spec_logmsg(logf, '[FOOOF] %s', strtrim(out)); end

if status ~= 0 || ~exist(outPath, 'file')
    error('spec_run_fooof_python:Fail', 'FOOOF python failed. See log.');
end

raw = fileread(outPath);
fooofOut = jsondecode(raw);
end