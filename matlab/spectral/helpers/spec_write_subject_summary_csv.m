function spec_write_subject_summary_csv(outPath, subjid, featGA, rclTable, logf)
% SPEC_WRITE_SUBJECT_SUMMARY_CSV  One-row subject-level summary file
% V 1.0.0
%
% Computes and writes TVI_alpha (temporal variability index) from the
% per-trial GA BI_pre sequence, plus optional GA circular-linear r_cl
% from the phase metric.
%
% TVI_alpha = MSSD / Var(BI_pre)
%   MSSD = mean square successive difference of BI_pre across trials
%   Range [0, 2]; near 0 = rigid/slowly-varying; near 2 = rapidly alternating.
%   Requires >= 3 valid trials.
%
% Inputs:
%   outPath  : full path to the output CSV (rows are appended, not overwritten)
%   subjid   : integer subject ID
%   featGA   : GA feature struct (must include bi_pre [nTrials x 1])
%   rclTable : table from spec_compute_phase_metric, or [] / empty table
%   logf     : (optional) MATLAB log file handle

if nargin < 5, logf = 1; end

% ---------------------------------------------------------------
% TVI_alpha from GA BI_pre sequence
% ---------------------------------------------------------------
if ~isfield(featGA, 'bi_pre')
    error('spec_write_subject_summary_csv:MissingBI', ...
        'featGA must contain a bi_pre field. Run spec_compute_interaction_metrics first.');
end

b      = featGA.bi_pre(:);    % [nTrials x 1]
b_ok   = b(~isnan(b));
K      = numel(b_ok);

tvi_alpha = nan;
mssd_bi   = nan;
var_bi    = nan;

if K >= 3
    diffs   = diff(b_ok);
    mssd_bi = mean(diffs .^ 2);
    var_bi  = var(b_ok, 0);    % unbiased (denominator K-1)

    if var_bi > eps
        tvi_alpha = mssd_bi / var_bi;
    else
        spec_logmsg(logf, '[SUBJ_SUM] sub-%03d: Var(BI_pre) ~ 0; TVI_alpha = NaN.', subjid);
    end
else
    spec_logmsg(logf, '[SUBJ_SUM][WARN] sub-%03d: only %d valid BI_pre trials (need >= 3).', subjid, K);
end

spec_logmsg(logf, '[SUBJ_SUM] sub-%03d  n=%d  TVI_alpha=%.4f  MSSD=%.4g  Var=%.4g', ...
    subjid, K, tvi_alpha, mssd_bi, var_bi);

% ---------------------------------------------------------------
% GA circular-linear r_cl (mean across channels)
% ---------------------------------------------------------------
ga_r_cl    = nan;
ga_p_perm  = nan;
n_chan_rcl = 0;

hasRcl = ~isempty(rclTable) && istable(rclTable) && ...
         ismember('r_cl', rclTable.Properties.VariableNames) && ...
         height(rclTable) > 0;

if hasRcl
    valid_r = rclTable.r_cl(~isnan(rclTable.r_cl));
    valid_p = rclTable.p_perm(~isnan(rclTable.p_perm));
    if ~isempty(valid_r)
        ga_r_cl    = mean(valid_r);
        ga_p_perm  = mean(valid_p);
        n_chan_rcl = numel(valid_r);
        spec_logmsg(logf, '[SUBJ_SUM] GA r_cl=%.4f  p_perm=%.4f  n_chans=%d', ...
            ga_r_cl, ga_p_perm, n_chan_rcl);
    end
end

% ---------------------------------------------------------------
% Build single-row table and write
% ---------------------------------------------------------------
T = table( ...
    int32(subjid),     int32(K),      tvi_alpha, ...
    mssd_bi,           var_bi,        ga_r_cl,   ...
    ga_p_perm,         int32(n_chan_rcl),         ...
    'VariableNames', { ...
        'subjid',          'n_trials_bi',   'tvi_alpha', ...
        'mssd_bi_pre',     'var_bi_pre',    'ga_phase_r_cl', ...
        'ga_phase_p_perm', 'n_chans_r_cl'  } ...
);

% Append rows so all subjects land in one CSV without overwriting
if exist(outPath, 'file')
    writetable(T, outPath, 'WriteMode', 'append', 'WriteVariableNames', false);
else
    writetable(T, outPath);
end

spec_logmsg(logf, '[SUBJ_SUM] Written -> %s', outPath);
end