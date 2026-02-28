function spec_write_chan_trial_csv(outPath, subjid, chanLabels, feat)
% Writes one row per (trial, chan)

nChan = numel(chanLabels);
[nChan, nTr] = size(feat.paf_cog_hz);
if nChan ~= nChan
    error('spec_write_chan_trial_csv:Shape', 'chan label count mitmatch');
end

trialCol = repelem((1:nTr)', nChan);
chanIdxCol = repmat((1:nChan)', nTr, 1);
chanLabCol = repmat(string(chanLabels(:)), nTr, 1);

T = table();
T.subjid = repmat(subjid, nChan * nTr, 1);
T.trial = trialCol;
T.chan_idx = chanIdxCol;
T.chan_label = chanLabCol;

% Helper to flatten [chan x trial] -> [nChan * nTr x 1] in trial-major order
flat = @(X) reshape(X, [nChan * nTr, 1]);

T.paf_cog_hz        = flat(feat.paf_cog_hz);
T.pow_slow_alpha    = flat(feat.pow_slow_alpha);
T.pow_fast_alpha    = flat(feat.pow_fast_alpha);
T.pow_alpha_total   = flat(feat.pow_alpha_total);
T.rel_slow_alpha    = flat(feat.rel_slow_alpha);
T.rel_fast_alpha    = flat(feat.rel_fast_alpha);
T.sf_ratio          = flat(feat.sf_ratio);
T.sf_logratio       = flat(feat.sf_logratio);
T.sf_balance        = flat(feat.sf_balance);
T.slow_alpha_frac   = flat(feat.slow_alpha_frac);

writetable(T, outPath);
end