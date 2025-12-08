function exp01_preproc(only_sub, doICA, doPost, interactiveClean)
% EXP01_PREPROC Wrapper for the 26ByBiosemi experiment (exp01)
%
% exp01_preproc(sonly_sub, doICA, doPost, interactiveClean)
%
% only_sub : 'all' deault, or specific ID (e.g., 'sub-01')
% doICA : true/false (default: false)
% doPost: true/false (default: false)
% interactiveClean : true/false (default: same as doICA)
%
% This calls the generic preproc_default() with experiment specific
% identifiers and configuration

% Set default values for optional parameters
    if nargin < 1 || isempty(only_sub), only_sub = 'all'; end
    if nargin < 2 || isempty(doICA), doICA = false; end
    if nargin < 3 || isempty(doPost), doPost = false; end
    if nargin < 4 || isempty(interactiveClean), interactiveClean = doICA; end

    % Experiment ID is used by config_paths (BIDS folder name)
    exp_id = '26ByBiosemi';
    cfg_file = 'exp01.json';

    preproc_default(exp_id, cfg_file, only_sub, doICA, doPost, interactiveClean);
end