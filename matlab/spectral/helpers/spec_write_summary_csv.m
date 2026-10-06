function spec_write_summary_csv(outPath, subjid, chanLabels, S)
% SPEC_WRITE_SUMMARY_CSV Subject x channel table (one row per channel)
% V 1.0.0
%
% Channel-space counterpart of the source pipeline's sub-XXX_source_ga.csv
% (subject x ROI): same column names, computed the same way (trial-mean
% powers / PSD / waveform, per-channel TVI, ITC, FOOOF).
%
% Inputs:
%   S : struct of [nChan x 1] numeric fields; other shapes are skipped.

if isstring(chanLabels), chanLabels = cellstr(chanLabels); end
nChan = numel(chanLabels);

T = table();
T.subjid = repmat(int32(subjid), nChan, 1);
T.chan_idx = (1:nChan)';
T.chan_label = string(chanLabels(:));

fns = sort(fieldnames(S));
for k = 1:numel(fns)
    v = S.(fns{k});
    if (isnumeric(v) || islogical(v) && numel(v) == nChan)
        T.(fns{k}) = double(v(:));
    end
end

writetable(T, outPath);
end