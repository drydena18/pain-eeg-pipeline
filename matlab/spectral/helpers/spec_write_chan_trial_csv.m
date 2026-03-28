function spec_write_chan_trial_csv(outPath, subjid, chanLabels, featChan)
% SPEC_WRITE_CHAN_TRIAL_CSV  One row per (trial, channel)
% V 1.1.0  -- bug fixes + new interaction metric and phase columns
%
% New columns vs V1.0:
%   bi_pre, lr_pre, cog_pre, psi_cog  (pre-stim sub-band metrics)
%   erd_slow, erd_fast, delta_erd     (ERD asymmetry family)
%   p5_flag                           (unstable-power trial flag)
%   phase_slow_rad                    (Hilbert instantaneous phase; NaN if unavailable)

% BUG FIX: original code assigned nChan = numel(chanLabels) then immediately
% overwrote it with [nChan, nTr] = size(...), making the guard always false.
nLabels = numel(chanLabels);                       % label count (for validation)
[nChan, nTr] = size(featChan.paf_cog_hz);

if nLabels ~= nChan
    error('spec_write_chan_trial_csv:Mismatch', ...  % BUG FIX: was 'mitmatch'
        'Channel label count (%d) != feature array rows (%d).', nLabels, nChan);
end

trialCol   = repelem((1:nTr)',    nChan);
chanIdxCol = repmat((1:nChan)',   nTr,   1);
chanLabCol = repmat(string(chanLabels(:)), nTr, 1);

% Helper: [chan x trial] -> [nChan*nTr x 1] in trial-major order
flat = @(X) reshape(X, [nChan * nTr, 1]);

T = table();
T.subjid      = repmat(subjid, nChan * nTr, 1);
T.trial       = trialCol;
T.chan_idx    = chanIdxCol;
T.chan_label  = chanLabCol;

% ------ Existing full-epoch spectral features ------
T.paf_cog_hz       = flat(featChan.paf_cog_hz);
T.pow_slow_alpha   = flat(featChan.pow_slow_alpha);
T.pow_fast_alpha   = flat(featChan.pow_fast_alpha);
T.pow_alpha_total  = flat(featChan.pow_alpha_total);
T.rel_slow_alpha   = flat(featChan.rel_slow_alpha);
T.rel_fast_alpha   = flat(featChan.rel_fast_alpha);
T.sf_ratio         = flat(featChan.sf_ratio);
T.sf_logratio      = flat(featChan.sf_logratio);
T.sf_balance       = flat(featChan.sf_balance);
T.slow_alpha_frac  = flat(featChan.slow_alpha_frac);

% ------ Pre-stimulus interaction metrics (new) ------
if isfield(featChan, 'bi_pre')
    T.bi_pre     = flat(featChan.bi_pre);
    T.lr_pre     = flat(featChan.lr_pre);
    T.cog_pre    = flat(featChan.cog_pre);
    T.psi_cog    = flat(featChan.psi_cog);
    T.erd_slow   = flat(featChan.erd_slow);
    T.erd_fast   = flat(featChan.erd_fast);
    T.delta_erd  = flat(featChan.delta_erd);
    T.p5_flag    = flat(featChan.p5_flag);
end

% ------ Slow-alpha Hilbert instantaneous phase (new; optional) ------
if isfield(featChan, 'phase_slow_rad')
    T.phase_slow_rad = flat(featChan.phase_slow_rad);
else
    T.phase_slow_rad = nan(nChan * nTr, 1);   % NaN sentinel: stage 09 not run
end

writetable(T, outPath);
end