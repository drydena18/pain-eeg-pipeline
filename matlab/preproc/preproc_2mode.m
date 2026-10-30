function out = preproc_2mode(exp_id, action, subjects_override)
% PREPROC_2MODE Two-mode preprocessing for ONE experiment.
% V 1.0.0
%
%   preproc_2mode('exp06')                      % unattended: Mode 1 -> Mode 2a
%                                               % (and Mode 2b for anyone already reviewed)
%   preproc_2mode('exp06', 'review')            % manual IC rejection, then Mode 2b
%   preproc_2mode('exp06', 'status')            % print the manifest only
%   preproc_2mode('exp06', 'run', [1 4 17])     % subset of subjects
%   preproc_2mode('exp06', 'rerun_stale')       % redo ICA wher Mode 1 config changed
%
% Call chain:
%   preproc_2mode -> pp_from_experiment (+ pp_hooks) -> pp_run / pp_review
%
% Same config (config/<exp>.json), same config_pahts, same per-subject
% folders and file names as expXX_preproc -> preproc_core. See
% two_mode/README.md for what differs

if nargin < 2 || isempty(action), action = "run"; end
if nargin < 3, subjects_override = []; end

here = fileparts(mfilename("fullpath"));
matlabDir = fileparts(here);
addpath(fullfile(here, 'two_mode'), fullfile(here, 'two_mode', 'stages'), fullfile(matlabDir, 'helpers'), fullfile(matlabDir, 'utils'));

% ----------------------------------
% EEGLAB init (EDIT ON EACH MACHINE)
% ----------------------------------
try
    if exist('eeglab', 'file') ~= 2
        addpath('/home/UWO/darsenea/Documents/matlab-toolboxes/eeglab2025.1.0');
    end
    if exist('eeglab', 'file') ~= 2
        error('EEGLAB not on path. Add EEGLAB folder to MATLAB path.');
    end
    if exist('pop_loadset', 'file') ~= 2 || exist('iclabel', 'file') ~= 2
        eeglab nogui;
    end
catch ME
    error('preproc_2mode:EEGLABInitFail', 'EEGLAB init failed: %s', ME.message);
end

pp = pp_from_experiment(exp_id, subjects_override);

switch lower(char(action))
    case 'run'
        out = pp_run(pp);
    case 'review'
        out = pp_manifest(pp, pp.subjects);
    case 'status'
        pp_ensure_dirs(pp);
        out = pp_manifest(pp, pp.subjects);
        disp(out);
    case 'rerun_stale'
        out = pp_run(pp, 'rerun_stale', true);
    otherwise
        error('preproc_2mode:Action', 'Unknown action "%s" (run | review | status | rerun_stale).', action);
end
end