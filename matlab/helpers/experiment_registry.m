function R = experiment_registry()
% Define all experiments here.
% raw_dirname: folder under RAW_ROOT that contains raw data
% out_dirname: folder name under PROJ_ROOT for outputs
% out_prefix:  filename prefix for outputs (optional)

R = struct();

% exp01 - 26ByBiosemi
R.exp01 = struct( ...
    'id', 'exp01', ...
    'raw_dirname', '26ByBiosemi', ...
    'out_dirname', '26ByBiosemi', ...
    'out_prefix', '26BB_62_' ...
    );

end