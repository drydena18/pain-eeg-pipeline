function exp07_preproc(subjects_override)
% EXP07_PREPROC Entrypoint for preprocessing experiment 07.
% V 1.0.0
%
% Call chain:
%   exp07_preproc.m -> preproc_default.m -> preproc_core.m

exp_id = "exp07";

% Paths
P = config_paths(exp_id);

% ----------------------------------
% EEGLAB init (EDIT ON EACH MACHINE)
% ----------------------------------
% If EEGLAB is already on the path, this will just no-op
try
    if exist('eeglab', 'file') ~= 2
        % EDIT if needed
        addpath('/home/UWO/darsenea/Documents/matlab-toolboxes/eeglab2025.1.0');
    end
    if exist('eeglab', 'file') ~= 2
        error('EEGLAB not on path. Add EEGLAB folder to MATLAB path.');
    end
    eeglab nogui;
catch ME
    error('exp07_preproc:EEGLABInitFail', 'EEGLAB init failed: %s', ME.message);
end

% -------------------------
% Load experiment JSON cfg
% -------------------------
cfgFileDir = fullfile('/home/UWO/darsenea/Documents/GitHub/pain-alpha-dynamics/config/');
cfg_path = fullfile(cfgFileDir, sprintf('%s.json', exp_id));

if ~exist(cfg_path, 'file')
    error('exp07_preproc:MissingJSON', 'Config JSON not found: %s', cfg_path);
end

cfg = load_cfg(cfg_path);

% Sanity prints
fprintf('cfg class: %s\n', class(cfg));
fprintf('exp id: %s\n', string(cfg.exp.id));
fprintf('raw pattern: %s\n', string(cfg.exp.raw.pattern));

% ------------
% Build paths
% ------------
P = config_paths(exp_id, cfg);

% ----------------------------------
% Run default normalizer + pipeline
% ----------------------------------
preproc_default(P, cfg, subjects_override);

end