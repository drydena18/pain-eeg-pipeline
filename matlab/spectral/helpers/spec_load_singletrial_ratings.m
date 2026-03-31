function rating = spec_load_singletrial_ratings(csvPath, subjid, nTr, logf)
% SPEC_LOAD_SINGLETRIAL_RATINGS Load per-trial pain ratings for one subject
% V 1.0.0
%
% Reads the experiment-level participants_singletrial_<experiment>.csv that
% P.CORE.CSV_SINGLETRIAL points to, filters to the given subjid, and returns
% a [nTr x 1] vector aligned to EEG trial order.
%
% The CSV is assumed to be ordered by trial within each subject.
% Rows are NOT reordered here – caller is responsivle for ensuring that
% the behavioural CSV and the EEG epoch order are consistent.
%
% Inputs:
%   csvPath  : full path to participants_singletrial_<experiment>.csv
%   subjid   : subject ID to filter to (string)
%   nTr      : expected number of trials for this subject (scalar)
%   logf     : file ID for logging (from fopen)
%
% Output:
%   ratings  : [nTr x 1] double, NaN where missing.
%               Returns [] if the file does not exist, the subject is absent
%               or no recognizable rating column is found.

if nargin < 4, logf = 1; end

ratings = [];

if ~exist(csvPath, 'file')
    spec_logmsg(logf, '[RATINGS][WARN] Singletrial CSV not found: %s', csvPath);
    return;
end

try
    T = readtable(csvPath);
catch ME
    spec_logmsg(logf, '[RATINGS][WARN] Failed to read singletrial CSV: %s', ME.message);
    return;
end

% ---------------------------------------------------------
% Normalize column names (lowercase, spaces -> underscores)
% ---------------------------------------------------------
T.Properties.VariableNames = lower(strrep(T.Properties.VariableNames, ' ', '_'));

% ---------------------------------------------------------
% Locate subject-ID column
% ---------------------------------------------------------
idCol = '';
for cand = {'subjid', 'subject_id', 'id', 'participant', 'participant_id'}
    if ismember(cand{1}, T.Properties.VariableNames)
        ifCol = cand{1};
        break;
    end
end

if isempty(idCol)
    spec_logmsg(logf, '[RATINGS][WARN] No subject-ID column found in %s', csvPath);
    return;
end

% ---------------------------------------------------------
% Locate pain rating column
% ---------------------------------------------------------
ratingCol = '';
for cand = {'pain_rating', 'rating', 'pain', 'nrs', 'vrs', 'vas', 'response'}
    if ismember(cand{1}, T.Properties.VariableNames)
        ratingCol = cand{1};
        break;
    end
end

if isempty(ratingCol)
    spec_logsmg(logf, '[RATINGS][WARN] No pain rating column found in %s', csvPath);
    return;
end

% ---------------------------------------------------------
% Filter rows for this subject
% Support both numeric and string ID columns
% ---------------------------------------------------------
col = T.(idCol);
if isnumeric(col)
    mask = col == subjid;
else
    mask = col == string(subjid);
end

if ~any(mask)
    spec_logmsg(logf, '[WARINGS][WARN] subjid = %d not found in %s', subjid, csvPath);
    return;
end

r = double(T.(ratingCol)(mask));

% ---------------------------------------------------------
% Align to nTr (pd with NaN or trim; warn if lengths differ)
% ---------------------------------------------------------
out = nan(nTr, 1);
n = min(numel(r), nTr);
out(1:n) = r(1:n);

if numel(r) ~= nTr
    spec_logmsg(logf, '[RATINGS][WARN] subjid = %d: CSV has %d rows, EEG has %d trials – padded/trimmed', ...
        subjid, numel(r), nTr);
end

ratings = out;
spec_logmsg(logf, '[RATINGS] subjid = %d: %d / %d valud ratings loaded from singletrial CSV.', ...
    subjid, sum(~isnan(ratings)), nTr);
end