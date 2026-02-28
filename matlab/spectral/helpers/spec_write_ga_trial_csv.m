function spec_write_ga_trial_csv(outPath, subjid, featGA, fooofOut)
% SPEC_WRITE_GA_TRIAL_CSV
% One row per trial (GA spectral features + FOOOF outputs if available)

nTr = numel(featGA.paf_cog_hz);

T = table();
T.subjid = repmat(subjid, nTr, 1);
T.trial  = (1:nTr)';

% --- GA spectral features ---
T.paf_cog_hz      = featGA.paf_cog_hz(:);
T.pow_slow_alpha  = featGA.pow_slow_alpha(:);
T.pow_fast_alpha  = featGA.pow_fast_alpha(:);
T.pow_alpha_total = featGA.pow_alpha_total(:);
T.rel_slow_alpha  = featGA.rel_slow_alpha(:);
T.rel_fast_alpha  = featGA.rel_fast_alpha(:);
T.sf_ratio        = featGA.sf_ratio(:);
T.sf_logratio     = featGA.sf_logratio(:);
T.sf_balance      = featGA.sf_balance(:);
T.slow_alpha_frac = featGA.slow_alpha_frac(:);

% --- FOOOF outputs (raw + filled + provenance) ---
hasFooof = (nargin >= 4) && ~isempty(fooofOut) && isfield(fooofOut, 'trials') && ~isempty(fooofOut.trials);

if hasFooof
    nAvail = min(nTr, numel(fooofOut.trials));

    % Always include these (raw)
    off = nan(nTr,1); ex = nan(nTr,1); r2 = nan(nTr,1); er = nan(nTr,1);
    acf = nan(nTr,1); apw = nan(nTr,1); abw = nan(nTr,1);

    % Optional "filled" + provenance
    acf_f = nan(nTr,1); apw_f = nan(nTr,1); abw_f = nan(nTr,1);
    acf_src = strings(nTr,1); apw_src = strings(nTr,1); abw_src = strings(nTr,1);
    alpha_found = nan(nTr,1);

    for t = 1:nAvail
        s = fooofOut.trials(t);

        if isfield(s,'aperiodic_offset'),   off(t) = s.aperiodic_offset; end
        if isfield(s,'aperiodic_exponent'), ex(t)  = s.aperiodic_exponent; end
        if isfield(s,'r2'),                r2(t)  = s.r2; end
        if isfield(s,'error'),             er(t)  = s.error; end

        % Raw alpha peak outputs from python (your fooof_bridge.py)
        if isfield(s,'alpha_cf'), acf(t) = s.alpha_cf; end
        if isfield(s,'alpha_pw'), apw(t) = s.alpha_pw; end
        if isfield(s,'alpha_bw'), abw(t) = s.alpha_bw; end

        % Filled outputs from spec_fill_fooof_alpha (if you ran it)
        if isfield(s,'alpha_cf_filled'), acf_f(t) = s.alpha_cf_filled; end
        if isfield(s,'alpha_pw_filled'), apw_f(t) = s.alpha_pw_filled; end
        if isfield(s,'alpha_bw_filled'), abw_f(t) = s.alpha_bw_filled; end

        if isfield(s,'alpha_cf_source'), acf_src(t) = string(s.alpha_cf_source); end
        if isfield(s,'alpha_pw_source'), apw_src(t) = string(s.alpha_pw_source); end
        if isfield(s,'alpha_bw_source'), abw_src(t) = string(s.alpha_bw_source); end

        if isfield(s,'alpha_found'), alpha_found(t) = s.alpha_found; end
    end

    T.fooof_offset   = off;
    T.fooof_exponent = ex;
    T.fooof_r2       = r2;
    T.fooof_error    = er;

    T.fooof_alpha_cf = acf;
    T.fooof_alpha_pw = apw;
    T.fooof_alpha_bw = abw;

    % Only add filled/provenance columns if they exist in output
    if any(~isnan(acf_f)) || any(acf_src ~= "")
        T.fooof_alpha_cf_filled = acf_f;
        T.fooof_alpha_cf_source = acf_src;
    end
    if any(~isnan(apw_f)) || any(apw_src ~= "")
        T.fooof_alpha_pw_filled = apw_f;
        T.fooof_alpha_pw_source = apw_src;
    end
    if any(~isnan(abw_f)) || any(abw_src ~= "")
        T.fooof_alpha_bw_filled = abw_f;
        T.fooof_alpha_bw_source = abw_src;
    end
    if any(~isnan(alpha_found))
        T.fooof_alpha_found = alpha_found;
    end
end

writetable(T, outPath);
end