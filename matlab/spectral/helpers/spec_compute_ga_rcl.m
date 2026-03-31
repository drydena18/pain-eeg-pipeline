function gaRcl = spec_compute_ga_rcl(rclTable, logf)
% SPEC_COMPUTE_GA_RCL  Grand-average circular-linear correlation across channels
% V 1.0.0
%
% Reduces the per-channel r_cl / p_perm table produced by
% spec_compute_phase_metric into a single subject-level summary by taking
% the arithmetic mean across channels that have valid (non-NaN) values.
%
% This is a simple aggregation step — no new statistics are computed here.
% The per-channel r_cl values themselves are computed in
% spec_compute_phase_metric via the Harrison-Kanji formula.
%
% Inputs:
%   rclTable : table from spec_compute_phase_metric, with columns
%                chan_idx, chan_label, r_cl, p_perm
%              Pass [] or an empty table if phase was not computed.
%   logf     : (optional) MATLAB log file handle
%
% Output:
%   gaRcl : struct
%             .ga_phase_r_cl    scalar  mean r_cl across valid channels (NaN if unavailable)
%             .ga_phase_p_perm  scalar  mean p_perm across valid channels (NaN if unavailable)
%             .n_chans_r_cl     integer number of channels with valid r_cl

if nargin < 2, logf = 1; end

gaRcl = struct( ...
    'ga_phase_r_cl',   nan, ...
    'ga_phase_p_perm', nan, ...
    'n_chans_r_cl',    0   ...
);

hasRcl = ~isempty(rclTable) && istable(rclTable) && ...
         ismember('r_cl',   rclTable.Properties.VariableNames) && ...
         ismember('p_perm', rclTable.Properties.VariableNames) && ...
         height(rclTable) > 0;

if ~hasRcl
    spec_logmsg(logf, '[GA_RCL] No rclTable available; GA r_cl = NaN.');
    return;
end

valid_r = rclTable.r_cl(~isnan(rclTable.r_cl));
valid_p = rclTable.p_perm(~isnan(rclTable.p_perm));

if isempty(valid_r)
    spec_logmsg(logf, '[GA_RCL] All per-channel r_cl values are NaN.');
    return;
end

gaRcl.ga_phase_r_cl   = mean(valid_r);
gaRcl.ga_phase_p_perm = mean(valid_p);
gaRcl.n_chans_r_cl    = numel(valid_r);

spec_logmsg(logf, '[GA_RCL] GA r_cl=%.4f  GA p_perm=%.4f  n_chans=%d', ...
    gaRcl.ga_phase_r_cl, gaRcl.ga_phase_p_perm, gaRcl.n_chans_r_cl);
end