function EEG = apply_bids_channel_labels(EEG, tsvPath, logf)
% APPLY_BIDS_CHANNEL_LABELS  Relabel EEG.chanlocs from a BIDS channels.tsv
% V 1.0.0
%
% Reads the 'name' column from a BIDS-format channels.tsv and applies those
% labels to EEG.chanlocs in row order. This is run before the .elp coordinate
% lookup so that label matching succeeds.
%
% The TSV must have a 'name' column (BIDS required field). Row count must
% match EEG.nbchan exactly; if it does not, the function warns and returns
% EEG unchanged so downstream stages are not silently corrupted.
%
% Inputs:
%   EEG     : EEGLAB EEG struct
%   tsvPath : char/string path to channels.tsv
%   logf    : (optional) log file handle (default: stdout)
%
% Output:
%   EEG     : EEG struct with updated EEG.chanlocs labels

if nargin < 3, logf = 1; end

tsvPath = char(string(tsvPath));

if ~exist(tsvPath, 'file')
    logmsg(logf, '[CHANTSV][WARN] File not found, skipping: %s', tsvPath);
    return;
end

% Read TSV — MATLAB's readtable handles tab delimiters and headers cleanly
try
    T = readtable(tsvPath, 'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
catch ME
    logmsg(logf, '[CHANTSV][WARN] Could not read channels.tsv: %s', ME.message);
    return;
end

% Normalise column names to lowercase for robustness
T.Properties.VariableNames = lower(T.Properties.VariableNames);

if ~ismember('name', T.Properties.VariableNames)
    logmsg(logf, '[CHANTSV][WARN] channels.tsv has no ''name'' column. Columns found: %s', ...
        strjoin(T.Properties.VariableNames, ', '));
    return;
end

nTsv = height(T);
nEEG = EEG.nbchan;

if nTsv ~= nEEG
    logmsg(logf, '[CHANTSV][WARN] Row count mismatch: TSV has %d rows, EEG has %d channels. Skipping relabel.', ...
        nTsv, nEEG);
    return;
end

labels = string(T.name);

% Check for duplicates before applying — duplicate labels would corrupt
% any downstream label-based lookups
if numel(unique(labels)) ~= numel(labels)
    logmsg(logf, '[CHANTSV][WARN] channels.tsv contains duplicate names. Skipping relabel.');
    return;
end

for ch = 1:nEEG
    EEG.chanlocs(ch).labels = char(labels(ch));
end

logmsg(logf, '[CHANTSV] Applied %d channel labels from: %s', nEEG, tsvPath);
end