function spec_write_ga_trial_csv(outPath, subjid, featGA, fooofOut)
% SPEC_WRITE_GA_TRIAL_CSV  One row per trial (GA spectral features + FOOOF)
% V 1.1.0  -- adds pre-stim interaction metrics and phase column
%
% New columns vs V1.0:
%   bi_pre, lr_pre, cog_pre, psi_cog
%   erd_slow, erd_fast, delta_erd, p5_flag
%   phase_slow_rad   (circular GA; NaN if stage 09 not run)

nTr = numel(featGA.paf_cog_hz);

T          = table();
T.subjid   = repmat(subjid, nTr, 1);
T.trial    = (1:nTr)';

% ------ Existing full-epoch features ------
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

% ------ Pre-stimulus interaction metrics (new) ------
if isfield(featGA, 'bi_pre')
    T.bi_pre    = featGA.bi_pre(:);
    T.lr_pre    = featGA.lr_pre(:);
    T.cog_pre   = featGA.cog_pre(:);
    T.psi_cog   = featGA.psi_cog(:);
    T.erd_slow  = featGA.erd_slow(:);
    T.erd_fast  = featGA.erd_fast(:);
    T.delta_erd = featGA.delta_erd(:);
    T.p5_flag   = featGA.p5_flag(:);
end

% ------ Slow-alpha Hilbert phase GA (circular mean across channels) ------
if isfield(featGA, 'phase_slow_rad')
    T.phase_slow_rad = featGA.phase_slow_rad(:);
else
    T.phase_slow_rad = nan(nTr, 1);
end

% ------ FOOOF outputs (raw + filled + provenance) ------
hasFooof = (nargin >= 4) && ~isempty(fooofOut) && ...
           isfield(fooofOut, 'trials') && ~isempty(fooofOut.trials);

if hasFooof
    nAvail = min(nTr, numel(fooofOut.trials));

    off   = nan(nTr,1); ex  = nan(nTr,1); r2  = nan(nTr,1); er = nan(nTr,1);
    acf   = nan(nTr,1); apw = nan(nTr,1); abw = nan(nTr,1);
    acf_f = nan(nTr,1); apw_f = nan(nTr,1); abw_f = nan(nTr,1);
    acf_src = strings(nTr,1); apw_src = strings(nTr,1); abw_src = strings(nTr,1);
    alpha_found = nan(nTr,1);

    for t = 1:nAvail
        s = fooofOut.trials(t);
        off(t)  = getScalarNum(s, 'aperiodic_offset');
        ex(t)   = getScalarNum(s, 'aperiodic_exponent');
        r2(t)   = getScalarNum(s, 'r2');
        er(t)   = getScalarNum(s, 'error');
        acf(t)  = getScalarNum(s, 'alpha_cf');
        apw(t)  = getScalarNum(s, 'alpha_pw');
        abw(t)  = getScalarNum(s, 'alpha_bw');
        acf_f(t) = getScalarNum(s, 'alpha_cf_filled');
        apw_f(t) = getScalarNum(s, 'alpha_pw_filled');
        abw_f(t) = getScalarNum(s, 'alpha_bw_filled');
        acf_src(t) = getScalarStr(s, 'alpha_cf_source');
        apw_src(t) = getScalarStr(s, 'alpha_pw_source');
        abw_src(t) = getScalarStr(s, 'alpha_bw_source');
        alpha_found(t) = getScalarNum(s, 'alpha_found');
    end

    T.fooof_offset   = off;
    T.fooof_exponent = ex;
    T.fooof_r2       = r2;
    T.fooof_error    = er;
    T.fooof_alpha_cf = acf;
    T.fooof_alpha_pw = apw;
    T.fooof_alpha_bw = abw;

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

% ---------------- Helpers ----------------
function v = getScalarNum(s, fn)
    if ~isstruct(s) || ~isfield(s, fn) || isempty(s.(fn)); v = nan; return; end
    raw = s.(fn);
    if iscell(raw)
        if isempty(raw), v = nan; return; end
        raw = raw{1};
    end
    if isstring(raw) || ischar(raw)
        tmp = str2double(raw); if isnan(tmp), v = nan; else, v = tmp; end; return;
    end
    if islogical(raw); v = double(raw); return; end
    if isnumeric(raw)
        raw = double(raw); if isempty(raw), v = nan; else, v = raw(1); end; return;
    end
    v = nan;
end

function s1 = getScalarStr(s, fn)
    if ~isstruct(s) || ~isfield(s, fn) || isempty(s.(fn)); s1 = ""; return; end
    try
        s1 = string(s.(fn));
        if numel(s1) > 1, s1 = s1(1); end
    catch
        s1 = "";
    end
end