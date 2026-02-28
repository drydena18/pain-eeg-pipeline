function exp01_spectral(subjects_override)
% EXP01_SPECTRAL Entrypoint for spectral feature extraction (Experiment 1).
% V 1.0.1
%
% Call chain:
%   exp01_spectral.m -> spectral_default.m -> spectral_core.m
%
% Usage:
%   exp01_spectral();           % uses cfg.exp.subjects (or participants.tsv fallback)
%   exp01_spectral([1 2 3 10]); % override subjects list

% ----------------------------
% Add project paths
% ----------------------------
thisFile = mfilename('fullpath');

% If file lives at: <repo>/matlab/spectral/exp01_spectral.m
% then repo root is 2 levels up from /matlab/spectral/
repoRoot = fileparts(fileparts(fileparts(thisFile)));

addpath(genpath(fullfile(repoRoot, 'matlab')));

% Prefer *only* spectral helpers to avoid name collisions with preproc helpers
specHelpers = fullfile(repoRoot, 'matlab', 'spectral', 'helpers');
if exist(specHelpers, 'dir')
    addpath(genpath(specHelpers));
end

% -------------------
% Inputs / defaults
% -------------------
if nargin < 1
    subjects_override = [];
end

exp_id = "exp01";

% ------------------------------------
% EEGLAB init (edit on each machine)
% ------------------------------------
try
    if exist('eeglab', 'file') ~= 2
        addpath('/home/UWO/darsenea/Documents/matlab-toolboxes/eeglab2025.1.0');
    end
    if exist('eeglab', 'file') ~= 2
        error('EEGLAB not on path. Add EEGLAB folder to MATLAB path.');
    end
    eeglab nogui;
catch ME
    error('exp01_spectral:EEGLABInitFail', 'EEGLAB init failed: %s', ME.message);
end

% --------------------------
% Load experiment JSON cfg
% --------------------------
cfgFileDir = fullfile(repoRoot, 'config');   % <repo>/config
cfg_path   = fullfile(cfgFileDir, sprintf('%s.json', exp_id));

if ~exist(cfg_path, 'file')
    error('exp01_spectral:MissingJSON', 'Config JSON not found: %s', cfg_path);
end

cfg = load_cfg(cfg_path);

% -----------------------------
% Build Paths (registry + cfg)
% -----------------------------
P = config_paths(exp_id, cfg);

% ------------------------------
% Run spectral defaults + core
% ------------------------------
spectral_default(P, cfg, subjects_override);

end