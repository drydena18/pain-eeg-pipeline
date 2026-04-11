function spec_write_chan_trial_csv(outPath, subjid, chanLabels, feat)
% SPEC_WRITE_CHAN_TRIAL_CSV Write one row per (trial, channel) to CSV
% V 2.0.0
%
% V2.0.0: dynamically writes every numeric [nChan x nTr] field present in
% feat, rather than a hardcoded V1-only column list. New fields (bi_pre,
% lr_pre, cog_pre, psi_cog, erd_slow, erd_fast, delta_erd, p5_flag,
% phase_slow_rad, and any future additions) are written automatically with
% no changes required here.
%
% Fields that are not [nChan x nTr] (wrong shape, non-numeric, empty) are
% silently skipped so the function is rubust to partially-computed runs.
%
% Column order: subjid, trial, chan_idx, chan_label, <all numeric fields
% sorted alphabetically for deterministic CSV structure across subjects>.
%
% Inputs:
%   outPath     : full CSV outputh path
%   subjid      : integer subject ID
%   chanLabels  : cell or string array of channel labels (length nChan)
%   feat        : struct of [nChan x nTr] feature arrays

if isstring(chanLabels), chanLabels = cellstr(chanLabels); end

% Determine nChan and nTr from the first valid [nChan x nTr] field
nChan = numel(chanLabels);
nTr = [];
fns = fieldnames(feat);

for k = 1:numel(fns)
    v = feat.(fns{k});
    if isnumeric(v) && ismatrix(v) && size(v, 1) == nChan && size(v, 2) > 0
        nTr = size(v, 2);
        break;
    end
end

if isempty(nTr)
    warning('spec_write_chan_trial_csv:NoValidField', ...
        'No valid [nChan x nTr] field found in feat for sub-%03d. CSV not written', subjid);
        return;
end

nRows = nChan * nTr;

% ---------------------------------------------------------------
% Index columns (always present)
% ---------------------------------------------------------------
T = table();
T.subjid = repmat(int32(subjid), nRows, 1);
T.trial = repelem((1:nTr)', nChan);
T.chan_idx = repmat((1:nChan)', nTr, 1);
T.chan_label = repmat(string(chanLabels(:)), nTr, 1);


% ---------------------------------------------------------------
% Dynamic feature columns
% Flatten each valid [nChan x nTr] field in column-major order so
% rows are grouped as (trial 1, all chans), (trial 2, all chans), ...
% which matches the trial/chan_idx index columns above
% ---------------------------------------------------------------
flat = @(X) reshape(X, [nRows, 1]);

% Sort field names for deterministic column ordering
sortedFns = sort(fns);

for k = 1:numel(sortedFns)
    fn = sortedFns{k};
    v = feat.(fn);

    % Skip non-numeric, wrong-shape, or empty fields
    if ~isnumeric(v) && ~islogical(v)
        continue;
    end
    if ~ismatrix(v) || size(v, 1) ~= nChan || size(v, 2) ~= nTr
        continue;
    end

    T.(fn) = flat(double(v));
end

writetable(T, outPath);
end