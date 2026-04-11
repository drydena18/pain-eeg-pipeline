function spec_write_ga_trial_csv(outPath, subjid, featGA, fooofOut)
% SPEC_WRITE_GA_TRIAL_CSV  Write one row per trial (GA spectral features + FOOOF)
% V 2.0.0
%
% V2.0.0: dynamically writes every numeric [nTr x 1] (or [1 x nTr]) field
% present in featGA, rather than a hardcoded V1-only column list.  New GA
% fields (bi_pre, lr_pre, cog_pre, psi_cog, erd_slow, erd_fast, delta_erd,
% p5_flag, phase_slow_rad, and any future additions) are written
% automatically with no changes required here.
%
% The FOOOF block (raw alpha peak + filled + provenance columns) is
% unchanged from V1.  It is appended after all GA fields.
%
% Column order: subjid, trial, <all numeric GA fields sorted
% alphabetically>, <FOOOF columns if available>.
%
% Inputs:
%   outPath  : full CSV output path
%   subjid   : integer subject ID
%   featGA   : struct of [nTr x 1] or [1 x nTr] feature arrays
%   fooofOut : (optional) struct from spec_run_fooof_python / spec_fill_fooof_alpha

if nargin < 4, fooofOut = struct(); end

% Determine nTr from the first valid field
fns = fieldnames(featGA);
nTr = [];
for k = 1:numel(fns)
    v = featGA.(fns{k});
    if isnumeric(v) && isvector(v) && numel(v) > 0
        nTr = numel(v);
        break;
    end
end

if isempty(nTr)
    warning('spec_write_ga_trial_csv:NoValidField', ...
        'No valid vector field found in featGA for sub-%03d. CSV not written.', subjid);
    return;
end

% ---------------------------------------------------------------
% Index columns
% ---------------------------------------------------------------
T = table();
T.subjid = repmat(int32(subjid), nTr, 1);
T.trial = (1:nTr)';

% ---------------------------------------------------------------
% Dynamic GA feature columns (sorted for deterministic ordering)
% ---------------------------------------------------------------
sortedFns = sort(fns);

for k = 1:numel(sortedFns)
    fn = sortedFns{k};
    v = featGA.(fn);

    if ~isnumeric(v) && ~islogical(v)
        continue;
    end

    if ~isvector(v) || numel(v) ~= nTr
        continue;
    end

    T.(fn) = double(v(:));
end

% ---------------------------------------------------------------
% FOOOF Outputs (raw + filled + provenance) - sunchanged from V1
% ---------------------------------------------------------------
hasFooof = ~isempty(fooofOut) && isstruct(fooofOut) && ...
    isfield(fooofOut, 'trials') && ~isempty(fooofOut.trials);

if hasFooof
    nAvail = min(nTr, numel(fooofOut.trials));

    off = nan(nTr,1); ex = nan(nTr,1); r2 = nan(nTr,1); er = nan(nTr,1);
    acf = nan(nTr,1); apw = nan(nTr,1); abw = nan(nTr,1);
    acf_f = nan(nTr,1); apw_f = nan(nTr,1); abw_f = nan(nTr,1);
    acf_src = strings(nTr,1); apw_src = strings(nTr,1); abw_src = strings(nTr,1);
    alpha_found = nan(nTr,1);

    for t = 1:nAvail
        s = fooofOut.trials(t);
        off(t) = getScalarNum(s, 'aperiodic_offset');
        ex(t)  = getScalarNum(s, 'aperiodic_exponent');
        r2(t)  = getScalarNum(s, 'r2');
        er(t)  = getScalarNum(s, 'error');
        acf(t) = getScalarNum(s, 'alpha_cf');
        apw(t) = getScalarNum(s, 'alpha_pw');
        abw(t) = getScalarNum(s, 'alpha_bw');
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

% ================================================================
% Local helpers (unchanged from V1)
% ================================================================
function v = getScalarNum(s, fn)
    if ~isstruct(s) || ~isfield(s, fn) || isempty(s.(fn))
        v = nan; return;
    end
    raw = s.(fn);
    if iscell(raw)
        if isempty(raw), v = nan; return; end
        raw = raw{1};
    end
    if isstring(raw) || ischar(raw)
        tmp = str2double(raw);
        if isnan(tmp), v = nan; else, v = tmp; end
        return;
    end
    if islogical(raw), v = double(raw); return; end
    if isnumeric(raw)
        raw = double(raw);
        if isempty(raw), v = nan; else, v = raw(1); end
        return;
    end
    v = nan;
end
 
function s1 = getScalarStr(s, fn)
    if ~isstruct(s) || ~isfield(s, fn) || isempty(s.(fn))
        s1 = ""; return;
    end
    try
        s1 = string(s.(fn));
        if numel(s1) > 1, s1 = s1(1); end
    catch
        s1 = "";
    end
end