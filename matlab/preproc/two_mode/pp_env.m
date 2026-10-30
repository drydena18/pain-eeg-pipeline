function env = pp_env(cfg)
%PP_ENV Software environment for provenance.
env.matlab = version;
env.eeglab = '';
try, env.eeglab = eeg_getversion; catch, end %#ok<NOCOM>
env.git = '';
if isfield(cfg.paths, 'code') && ~isempty(cfg.paths.code)
    [st, out] = system(sprintf('git -C "%s" rev-parse HEAD', cfg.paths.code));
    if st == 0, env.git = strtrim(out); end
end
env.host = '';
try, env.host = char(java.net.InetAddress.getLocalHost.getHostName); catch, end %#ok<NOCOM>
end