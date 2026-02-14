function exp01_preproc(subjects_override)
% EXP01_PREPROC Entrypoint for preprocessing experiment 01.
<<<<<<< HEAD
% V 1.0.3
%
% Call chain:
%   exp01_preproc.m -> preproc_default.m -> preproc_core.m
%
% Usage:    
%   exp01_preproc();            -> runs cfg.exp.subjects
%   exp01_preproc([1 2 3 10]);  -> override subjects list

% Initialize working directory
thisFile = filename('fullpath');
projRoot = fileparts(fileparts(thisFile));
addpath(genpath(projRoot));

%addpath(genpath('/path/to/pipeline/code')); % contains preproc_helpers.m


if nargin < 1
    subjects_override = [];
end

exp_id = "exp01";
=======
% V 1.0.2
%
% Call chain:
%   exp01_preproc.m -> preproc_default.m -> preproc_core.m

exp_id = "exp01";

% Paths
P = config_paths(exp_id);
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910

% ----------------------------------
% EEGLAB init (EDIT ON EACH MACHINE)
% ----------------------------------
% If EEGLAB is already on the path, this will just no-op
try
    if exist('eeglab', 'file') ~= 2
        % EDIT if needed
        %addpath('/path/to/eeglab');
    end
    if exist('eeglab', 'file') ~= 2
        error('EEGLAB not on path. Add EEGLAB folder to MATLAB path.');
    end
    eeglab nogui;
catch ME
    error('exp01_preproc:EEGLABInitFail', 'EEGLAB init failed: %s', ME.message);
end

<<<<<<< HEAD
% -------------------------
% Load experiment JSON cfg
% -------------------------
cfgFileDir = fullfile('/Users/drydena18/Desktop/pain-eeg-pipeline/config');
cfg_path = fullfile(cfgFileDir, sprintf('%s.json', exp_id));

=======
% ------------------------
% Load experiment JSON cfg
% ------------------------
cfg_path = fullfile(P.RESOURCE, sprintf('%s.json', exp_id));
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910
if ~exist(cfg_path, 'file')
    error('exp01_preproc:MissingJSON', 'Config JSON not found: %s', cfg_path);
end

<<<<<<< HEAD
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
=======
cfg = jsondecode(fileread(cfg_path));
disp(cfg.exp.id);
disp(cfg.exp.raw.pattern);

if nargin >= 1 && ~isempty(subjects_override)
    cfg.exp.subjects = subjects_override(:);
end

% ---------------------------------
% Run default normalizer + pipeline
% ---------------------------------
preproc_default(exp_id, P, cfg);
>>>>>>> bc259a06836bf8e2df96e96bf5570bd527991910

end