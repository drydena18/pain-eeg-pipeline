function out = spec_run_fooof_python(f, Prows, F, fooEnv, outTmp, subjid, logf)
% SPEC_RUN_FOOOF_PYTHON  Fit FOOOF/specparam to each row via fooof_bridge.py
% V 2.0.0
%
% The bridge calls python/source/src_fooof.src_fit_fooof, the same function
% the source pipeline uses, with the shared features.fooof parameters.
%
% V2.0.0 changes vs V1.x:
%   - Rows are channels (subject x channel trial-mean whole-epoch PSD),
%     not trials of the channel-GA PSD.
%   - Parameters come from cfg.spectral.feat.fooof (shared with source);
%     only python_exe / script_path come from cfg.spectral.fooof.
%   - Plain numeric CSV in / out; returns a struct of column vectors.
%
% Inputs:
%   f      : [1 x nF] Hz
%   Prows  : [nRows x nF] PSD
%   F      : cfg.spectral.feat
%   fooEnv : cfg.spectral.fooof (python_exe, script_path)
%   outTmp : scratch folder
%
% Output:
%   out : struct with fooof_offset, fooof_exponent, fooof_knee, fooof_r2,
%         fooof_error, fooof_alpha_cf, fooof_alpha_pw, fooof_alpha_bw,
%         each [nRows x 1]

names = {'fooof_offset', 'fooof_exponent', 'fooof_knee', 'fooof_r2', 'fooof_error', ...
         'fooof_alpha_cf', 'fooof_alpha_pw', 'fooof_alpha_bw'};

nRows = size(Prows, 1);
if size(Prows, 2) ~= numel(f)
    error('spec_run_fooof_python:Shape', 'Prows must be [nRows x nFreq].');
end
spec_ensure_dir(outTmp);

freqPath = fullfile(outTmp, sprintf('sub-%03d_fooof_freqs.csv', subjid));
psdPath  = fullfile(outTmp, sprintf('sub-%03d_fooof_psd.csv', subjid));
cfgPath  = fullfile(outTmp, sprintf('sub-%03d_fooof_cfg.json', subjid));
outPath  = fullfile(outTmp, sprintf('sub-%03d_fooof_out.csv', subjid));

writematrix(double(f(:)'), freqPath);
writematrix(double(Prows), psdPath);

pcfg = struct();
pcfg.fooof = F.fooof;
pcfg.alpha_band_hz = F.bands.alpha_hz(:)';
fid = fopen(cfgPath, 'w');
fwrite(fid, jsonencode(pcfg));
fclose(fid);

pyexe = "python3";
if isfield(fooEnv, 'python_exe') && strlength(string(fooEnv.python_exe)) > 0
    pyexe = string(fooEnv.python_exe);
end
script = string(fooEnv.script_path);
if strlength(script) == 0
    error('spec_run_fooof_python:MissingScript', 'cfg.spectral.fooof.script_path is required.');
end

if exist(outPath, 'file'), delete(outPath); end
cmd = sprintf('"%s" "%s" --freq "%s" --psd "%s" --cfg "%s" --out "%s" 2>&1', ...
    pyexe, script, freqPath, psdPath, cfgPath, outPath);
spec_logmsg(logf, '[FOOOF] CMD: %s', cmd);
[status, txt] = system(cmd);
if ~isempty(txt), spec_logmsg(logf, '[FOOOF] %s', strtrim(txt)); end
if status ~= 0 || ~exist(outPath, 'file')
    error('spec_run_fooof_python:Fail', 'FOOOF bridge failed (status %d). See log.', status);
end

M = readmatrix(outPath, 'NumHeaderLines', 1);
if size(M, 1) ~= nRows || size(M, 2) ~= numel(names)
    error('spec_run_fooof_python:BadOutput', 'Bridge output is %s, expected [%d %d].', ...
        mat2str(size(M)), nRows, numel(names));
end
out = struct();
for k = 1:numel(names)
    out.(names{k}) = M(:, k);
end
end
