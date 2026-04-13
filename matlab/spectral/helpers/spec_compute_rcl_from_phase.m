function rclTable = spec_compute_rcl_from_phase(phaseMat, ratings, chanLabels, phCfg, logf, P, subjid)
% SPEC_COMPUTE_RCL_FROM_PHASE Per-channel circular-linear correlation
% V 1.0.0
%
% Computes the Harrison-Kanji circular linear correlation coefficient (r_cl)
% between slow_alpha phase at stimulus onst and pain ratings, together with 
% a permutation-based p-value, for every channel.
%
% This is a pure computational function: it has no EEG struct or filtering
% dependancy and can be called whether phaseMat was obtained from the
% pre-computed 09_hilbert stage or from the inline Hilbert path.
%
% Inputs:
%   phaseMat    : [nChan x nTrials] instantaneous phase at t = 0 (radians, -pi to pi)
%   ratings     : [nTrials x 1] pain ratings (may contain NaN)
%   chanLabels  : cell or string array of channel labels, length nChan
%   phCfg       : cfg.spectral.phase; reads .n_permuations (default 5000)
%   logf        : log file handle for spec_logmsg
%
% Output:
%   rclTable : table with columns:
%                   chan_idx    : integer channel index
%                   chan_label  : string channel label
%                   r_cl        : double circular-linear correlation [0, 1]
%                   p_perm      : double one-sided permutation p-value
%                 Returns an empty table() is ratings is [] or all-NaN, or
%                 if nChan == 0.
%
% Notes on r_cl:
%   The Harrison-Kanji formula (1999):
%       r_cl = sqrt( (r_cs^2 + r_rc^2 - 2*r_cs*r_cc*r_sc) / (1 - r_sc^2) )
%   where r_cs = corr(y, sin(phi)), r_cc = corr(y, cos(phi)),
%         r_sc = corr(sin(phi), cos(phi)).
%   r_cl is always non-negative; it measures association strength but 
%   not direction. Directionality can be inspected via the sign of r_cs.
%
%   Permuation p-value uses the Phipson & Smyth (2010) formula
%       p = (count + 1) / (nPerm + 1)
%   to guarantee p > 0.

if nargin < 5, logf = 1; end

rclTable = table();
try
    rt = readtable(P.CORE.CSV_SINGLETRIAL, ...
        'VariableNamingRule', 'preserve');
    % Normalize column names
    rt.Properties.VariableNames = regexprep(...
        rt.Properties.VariableNames, '[^a-zA-Z0-9_]', '_');

    % Discover the subject column by case-insensitive matching 
    % mirrors the logic in extract_subject_column so both stay in sync
    vn = lower(string(rt.Properties.VariableNames));
    candidates = ["subjid", "participant_id", "id", "subject"];
    subColName = "";
    for c = candidates
        hit = find(vn == c, 1);
        if ~isempty(hit)
            subColName = rt.Properties.VariableNames{hit};
            break;
        end
    end

    if strlength(subColName) == 0
        spec_logmsg(logf, '[PHASE][WARN] No subject column found in ratings CSV. Columns: %s',  ...
            strjoin(rt.Properties.VariableNames, ', '));
        error('no_subject_col');
    end

    % Filter to this subject. Handle numeric vs string storage
    colVals = rt.(subColName);
    if isnumeric(colVals)
        mask = colVals == subjid;
    else
        mask = str2double(regexprep(string(colVals), '[^0-9]', '')) == subjid;
    end

    rcl_table = rt(mask, :);

    if isempty(rcl_table)
        spec_logmsg(logf, '[PHASE][WARN] Ratings CSV loaded but no rows match subjid = %d', subjid);
    else
        spec_logmsg(logf, '[PHASE] Loaded %d rating rows for subjid = %d', height(rcl_table), subjid);
    end

catch ME
    spec_logmsg(logf, '[PHASE][WARN] Could not load ratings: %s', ME.message);
end

nChan = size(phaseMat, 1);

if nChan == 0
    return;
end

if isempty(ratings) || all(isnan(ratings))
    spec_logmsg(logf, '[RCL][WARN] No valid ratings provided; skipping r_cl computation.');
    return;
end

% --------------------------------------------------
% Number of permutations
% --------------------------------------------------
nPerm = 5000;
if isfield(phCfg, 'n_permutations') && ~isempty(phCfg.n_permutations)
    nPerm = phCfg.n_permutations;
end

ratings = ratings(:); % ensure column

r_cl_vec = nan(nChan, 1);
p_perm_vec = nan(nChan, 1);

for ch = 1:nChan
    phi = phaseMat(ch, :)'; % [nTr x 1]
    ok = ~isnan(phi) & ~isnan(ratings);

    if sum(ok) < 10
        spec_logmsg(logf, '[RCL][WARN] Chan %d: only %d valid trials (need >= 10); skipping.', ...
            ch, sum(ok));
            continue;
    end

    r_cl_vec(ch) = circ_lin_corr(phi(ok), ratings(ok));
    p_perm_vec(ch) = permutation_p(phi(ok), ratings(ok), nPerm);
end

rclTable = table( ...
    (1:nChan)', string(chanLabels(:)), r_cl_vec, p_perm_vec, ...
    'VariableNames', {'chan_idx', 'chan_label', 'r_cl', 'p_perm'});

nValid = sum(~isnan(r_cl_vec));
spec_logmsg(logf, '[RCL] DONE: %d channels with valid r_cl. Range [%.3f, %.3f], median p = %.3f.', ...
    nValid, ...
    min(r_cl_vec, [], 'omitnan'), ...
    max(r_cl_vec, [], 'omitnan'), ...
    median(p_perm_vec, 'omitnan'));
end

% ==========================================================
% Local: Harrison-Kanji circular-linear correlation (1999)
% ==========================================================
function r = circ_lin_corr(phi, y)
s = sin(phi);
c = cos(phi);
r_cs = corr(y, s, 'type', 'Pearson');
r_cc = corr(y, c, 'type', 'Pearson');
r_sc = corr(s, c, 'type', 'Pearson');
denom = 1 - r_sc^2;
if abs(denom) < eps
    r = 0; return;
end
r2 = (r_cs^2 + r_cc^2 - 2*r_cs*r_cc*r_sc) / denom;
r = sqrt(max(r2, 0));
end

% ==========================================================
% Local: permutation p-value (onde-sided upper tail, H0: r_cl = 0)
% Phipson & Smyth (2010) formula: p = (count + 1) / (nPerm + 1)
% ==========================================================
function p = permutation_p(phi, y, nPerm)
obs = circ_lin_corr(phi, y);
nT = numel(phi);
count = 0;
for k = 1:nPerm
    phi_s = phi(randperm(nT));
    count = count + (circ_lin_corr(phi_s, y) >= obs);
end
p = (count + 1) / (nPerm + 1);
end