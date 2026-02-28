function fooofOut = spec_fill_fooof_alpha(fooofOut, f, gaPxx, featGA, alphaBandHz)
% SPEC_FILL_FOOOF_ALPHA
% Adds filled alpha peak params when FOOOF returned NaNs.
%
% Inputs:
%   fooofOut     : struct from jsondecode(fooof_out.json), with fooofOut.trials
%   f            : [1 x nFreq]
%   gaPxx        : [nFreq x nTrials] GA PSD per trial (power)
%   featGA       : GA spectral features (used for optional fallback)
%   alphaBandHz  : [8 12] typically
%
% Output:
%   fooofOut.trials(t) gets:
%       alpha_cf_filled, alpha_pw_filled, alpha_bw_filled
%       alpha_cf_source, alpha_pw_source, alpha_bw_source
%       alpha_found

if nargin < 5 || isempty(alphaBandHz)
    alphaBandHz = [8 12];
end

if isempty(fooofOut) || ~isfield(fooofOut,'trials') || isempty(fooofOut.trials)
    return;
end

% Ensure shapes
if size(gaPxx,1) ~= numel(f)
    error('spec_fill_fooof_alpha:Shape', 'gaPxx must be [nFreq x nTrials] to match f.');
end

nTr = min(size(gaPxx,2), numel(fooofOut.trials));

idxA = (f >= alphaBandHz(1)) & (f <= alphaBandHz(2));
fA = f(idxA);

% ---- Pre-create the new fields for ALL trials (safe) ----
% Doing it this way avoids structure mismatch issues later.
for t = 1:numel(fooofOut.trials)
    if ~isfield(fooofOut.trials(t), 'alpha_cf_filled'), fooofOut.trials(t).alpha_cf_filled = nan; end
    if ~isfield(fooofOut.trials(t), 'alpha_pw_filled'), fooofOut.trials(t).alpha_pw_filled = nan; end
    if ~isfield(fooofOut.trials(t), 'alpha_bw_filled'), fooofOut.trials(t).alpha_bw_filled = nan; end

    if ~isfield(fooofOut.trials(t), 'alpha_cf_source'), fooofOut.trials(t).alpha_cf_source = ""; end
    if ~isfield(fooofOut.trials(t), 'alpha_pw_source'), fooofOut.trials(t).alpha_pw_source = ""; end
    if ~isfield(fooofOut.trials(t), 'alpha_bw_source'), fooofOut.trials(t).alpha_bw_source = ""; end

    if ~isfield(fooofOut.trials(t), 'alpha_found'), fooofOut.trials(t).alpha_found = 0; end
end

% ---- Fill trial-by-trial ----
for t = 1:nTr
    s = fooofOut.trials(t);

    % Read raw FOOOF alpha if present
    raw_cf = getFieldOrNaN(s, 'alpha_cf');
    raw_pw = getFieldOrNaN(s, 'alpha_pw');
    raw_bw = getFieldOrNaN(s, 'alpha_bw');

    % If raw exists, mark as filled=raw and move on
    if ~isnan(raw_cf)
        fooofOut.trials(t).alpha_cf_filled = raw_cf;
        fooofOut.trials(t).alpha_cf_source = "fooof";
        fooofOut.trials(t).alpha_found = 1;
    end
    if ~isnan(raw_pw)
        fooofOut.trials(t).alpha_pw_filled = raw_pw;
        fooofOut.trials(t).alpha_pw_source = "fooof";
    end
    if ~isnan(raw_bw)
        fooofOut.trials(t).alpha_bw_filled = raw_bw;
        fooofOut.trials(t).alpha_bw_source = "fooof";
    end

    % If CF still missing, do fallback from GA PSD in alpha band
    if isnan(fooofOut.trials(t).alpha_cf_filled)
        y = gaPxx(idxA, t);

        if all(isnan(y)) || numel(y) < 3
            % last resort: use PAF CoG (still gives you “where alpha is”)
            if isfield(featGA,'paf_cog_hz') && numel(featGA.paf_cog_hz) >= t
                fooofOut.trials(t).alpha_cf_filled = featGA.paf_cog_hz(t);
                fooofOut.trials(t).alpha_cf_source = "paf_cog";
                fooofOut.trials(t).alpha_found = 1;
            end
        else
            % simple peak-pick within alpha band
            [pk, idx] = max(y);
            cf = fA(idx);

            % peak “power” (not dB) — this is NOT the same units as FOOOF peak PW,
            % but it’s a useful QC fallback.
            pw = pk;

            % crude bandwidth: width at half-max (Hz)
            half = pk / 2;
            left = find(y(1:idx) <= half, 1, 'last');
            right = idx - 1 + find(y(idx:end) <= half, 1, 'first');
            if isempty(left), left = 1; end
            if isempty(right), right = numel(y); end
            bw = fA(right) - fA(left);

            fooofOut.trials(t).alpha_cf_filled = cf;
            fooofOut.trials(t).alpha_pw_filled = pw;
            fooofOut.trials(t).alpha_bw_filled = bw;

            fooofOut.trials(t).alpha_cf_source = "psd_peak";
            fooofOut.trials(t).alpha_pw_source = "psd_peak";
            fooofOut.trials(t).alpha_bw_source = "fwhm";
            fooofOut.trials(t).alpha_found = 1;
        end
    end
end

end

function v = getFieldOrNaN(s, fn)
if isfield(s, fn)
    v = double(s.(fn));
else
    v = nan;
end
end