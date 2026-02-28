function spec_write_ga_trial_csv(outPath, subjid, featGA, fooofOut)
% One row per trial (GA features + FOOOF outputs if available)

nTr = numel(featGA.paf_cog_hz);
T = table();
T.subjid = repmat(subjid, nTr, 1);
T.trial = (1:nTr)';

T.paf_cog_hz        = featGA.paf_cog_hz;
T.pow_slow_alpha    = featGA.pow_slow_alpha;
T.pow_fast_alpha    = featGA.pow_fast_alpha;
T.pow_alpha_total   = featGA.pow_alpha_total;
T.rel_slow_alpha    = featGA.rel_slow_alpha;
T.rel_fast_alpha    = featGA.rel_fast_alpha;
T.sf_ratio          = featGA.sf_ratio;
T.sf_logratio       = featGA.sf_logratio;
T.sf_balance        = featGA.sf_balance;
T.slow_fast_frac    = featGA.slow_fast_frac;

% FOOOF outputs
if nargin >= 3 && ~isempty(fooofOut) && ifield(fooofOut, 'trials') && ~isempty(fooofOut.trials)
    % Map trial -> fields; assume fooofOut.trials is struct array length nTr
    ex = nan(nTr, 1); off = nan(nTr, 1); r2 = nan(nTr, 1); err = nan(nTr, 1);
    acf = nan(nTr, 1); apw = nan(nTr, 1); abw = nan(nTr, 1);

    for t = 1:min(nTr, numel(fooofOut.trials))
        s = fooofOut.trials(t);
        if isfield(s, 'aperiodic_exponent'), ex(t) = s.aperiodic_exponent; end
        if isfield(s, 'aperiodic_offset'), off(t) = s.aperiodic_offset; end
        if isfield(s, 'r2'), r2(t) = s.r2; end
        if isfield(s, 'error'), err(t) = s.error; end

        if isfield(s, 'alpha_peak_cf'), acf(t) = s.alpha_peak_cf; end
        if isfield(s, 'alpha_peak_pw'), apw(t) = s.alpha_peak_pw; end
        if isfield(s, 'alpha_peak_bw'), abw(t) = s.alpha_peak_bw; end
    end

    T.fooof_offset = off;
    T.fooof_exponent = ex;
    T.fooof_r2 = r2;
    T.fooof_error = err;
    T.fooof_alpha_peak_cf = acf;
    T.fooof_alpha_peak_pw = apw;
    T.fooof_alpha_peak_bw = abw;
end

writetable(T, outPath);
end