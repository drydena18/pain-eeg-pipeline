function ratings = spec_load_singletrial_ratings(csvPath, subjid, nTr, logf)
% SPEC_LOAD_SINGLETRIAL_RATINGS  Load per-trial pain ratings for one subject
% V 1.1.0
%
% Reads the experiment-level participants_singletrial_<experiment>.csv that
% P.CORE.CSV_SINGLETRIAL points to, filters to the given subjid, and returns
% a [nTr x 1] vector aligned to EEG trial order.
%
% V1.1.0 fix — column detection failure caused by readtable sanitisation:
%   readtable's default VariableNamingRule='modify' runs its own identifier
%   sanitisation (dots/spaces/carets -> underscores, leading digits prefixed,
%   etc.) BEFORE our code sees the names.  Applying strrep on top of that
%   produced double-garbled strings that never matched the candidate lists.
%
%   Fix: read with 'VariableNamingRule','preserve' (R2019b+) so MATLAB does
%   not touch the raw CSV header, then apply a single full normalisation:
%     1. strip leading/trailing whitespace
%     2. lowercase
%     3. replace any run of non-alphanumeric characters with one underscore
%     4. strip leading/trailing underscores
%
%   When no match is found the normalised column names are logged so the
%   problem is self-diagnosing without needing to open the CSV.
%
% Inputs:
%   csvPath : full path to participants_singletrial_<experiment>.csv
%             (i.e. P.CORE.CSV_SINGLETRIAL from config_paths.m)
%   subjid  : integer subject ID to filter on
%   nTr     : expected number of EEG trials — used for length validation
%   logf    : (optional) log file handle for spec_logmsg
%
% Output:
%   ratings : [nTr x 1] double, NaN where missing.
%             Returns [] if the file does not exist, the subject is absent,
%             or no recognisable rating column is found.

if nargin < 4, logf = 1; end

ratings = [];

if ~exist(csvPath, 'file')
    spec_logmsg(logf, '[RATINGS][WARN] Singletrial CSV not found: %s', csvPath);
    return;
end

% ---------------------------------------------------------------
% Read with preserved raw header (VariableNamingRule = 'preserve')
% Falls back to default read for MATLAB < R2019b.
% ---------------------------------------------------------------
try
    T = readtable(csvPath, 'TextType', 'string', 'VariableNamingRule', 'preserve');
catch
    try
        T = readtable(csvPath, 'TextType', 'string');
    catch ME
        spec_logmsg(logf, '[RATINGS][WARN] Failed to read singletrial CSV: %s', ME.message);
        return;
    end
end

% ---------------------------------------------------------------
% Full column-name normalisation applied to the preserved names.
% Examples:
%   "Subject ID"   -> "subject_id"
%   "pain.rating"  -> "pain_rating"
%   "NRS^2 Score"  -> "nrs_2_score"
%   "  VAS  "      -> "vas"
% ---------------------------------------------------------------
rawNames  = T.Properties.VariableNames;
normNames = normalize_colnames(rawNames);
T.Properties.VariableNames = normNames;

% ---------------------------------------------------------------
% Locate subject-ID column
% ---------------------------------------------------------------
idCands = {'subjid', 'subject_id', 'subjectid', 'id', ...
           'participant', 'participant_id', 'participantid', ...
           'sub', 'sub_id'};
idCol = find_column(normNames, idCands);

if isempty(idCol)
    spec_logmsg(logf, '[RATINGS][WARN] No subject-ID column found in %s', csvPath);
    spec_logmsg(logf, '[RATINGS][WARN]   Columns present (normalised): %s', ...
        strjoin(normNames, ', '));
    return;
end

% ---------------------------------------------------------------
% Locate pain-rating column
% ---------------------------------------------------------------
ratingCands = {'pain_rating', 'painrating', 'rating', 'pain', ...
               'nrs', 'nrs_score', 'vrs', 'vrs_score', ...
               'vas', 'vas_score', 'response', 'score'};
ratingCol = find_column(normNames, ratingCands);

if isempty(ratingCol)
    spec_logmsg(logf, '[RATINGS][WARN] No pain-rating column found in %s', csvPath);
    spec_logmsg(logf, '[RATINGS][WARN]   Columns present (normalised): %s', ...
        strjoin(normNames, ', '));
    return;
end

spec_logmsg(logf, '[RATINGS] id_col="%s"  rating_col="%s"', idCol, ratingCol);

% ---------------------------------------------------------------
% Filter rows for this subject
% Supports numeric and string ID columns.
% Also tries the "sub-NNN" format in case the CSV uses BIDS notation.
% ---------------------------------------------------------------
col = T.(idCol);
if isnumeric(col) || isinteger(col)
    mask = col == subjid;
else
    mask = (col == string(subjid)) | ...
           (col == sprintf('sub-%03d', subjid));
end

if ~any(mask)
    spec_logmsg(logf, '[RATINGS][WARN] subjid=%d not found in column "%s" of %s', ...
        subjid, idCol, csvPath);
    return;
end

r = double(T.(ratingCol)(mask));

% ---------------------------------------------------------------
% Align to nTr (pad with NaN or trim with warning)
% ---------------------------------------------------------------
out = nan(nTr, 1);
n   = min(numel(r), nTr);
out(1:n) = r(1:n);

if numel(r) ~= nTr
    spec_logmsg(logf, ...
        '[RATINGS][WARN] subjid=%d: CSV has %d rows, EEG has %d trials — padded/trimmed.', ...
        subjid, numel(r), nTr);
end

ratings = out;
spec_logmsg(logf, '[RATINGS] subjid=%d: %d / %d valid ratings loaded.', ...
    subjid, sum(~isnan(ratings)), nTr);
end

% ================================================================
%  Local: normalise cell array of column names to safe lowercase
%  identifiers with underscore separators.
% ================================================================
function out = normalize_colnames(names)
out = cell(size(names));
for k = 1:numel(names)
    s = strtrim(names{k});           % strip surrounding whitespace
    s = lower(s);                    % lowercase
    s = regexprep(s, '[^a-z0-9]+', '_');  % non-alphanumeric runs -> _
    s = regexprep(s, '^_+|_+$', '');      % strip leading/trailing _
    if isempty(s)
        s = sprintf('col%d', k);
    end
    out{k} = s;
end
end

% ================================================================
%  Local: return the first candidate present in normNames, else ''.
% ================================================================
function col = find_column(normNames, candidates)
col = '';
for k = 1:numel(candidates)
    if ismember(candidates{k}, normNames)
        col = candidates{k};
        return;
    end
end
end