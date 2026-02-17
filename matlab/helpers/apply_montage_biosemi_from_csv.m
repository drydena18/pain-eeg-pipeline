function EEG = apply_montage_biosemi_from_csv(P, cfg, EEG, logf, subjid, LOGS)
% Expects csv columns: raw_label, std_label
% Writes audit file: LOGS/sub-###_channelmap_applied.tsv

% Resolve csv path
csvPath = "";
if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'csv')
    csvPath = strtrim(string(cfg.exp.montage.csv));
end
if strlength(csvPath) == 0
    error('Montage enabled but cfg.exp.montage.csv is missing.');
end

% Prefer resource folder copy if relative
csvPath2 = fullfile(string(P.RESOURCE), char(csvPath));
if isfile(csvPath2)
    csvPath = string(csvPath2);
end

logmsg(logf, '[MONTAGE][DEBUG] P.RESOURCE = "%s"', string(P.RESOURCE));
logmsg(logf, '[MONTAGE][DEBUG] cfg.exp.montage.csv = "%s"', string(cfg.exp.montage.csv));
logmsg(logf, '[MONTAGE][DEBUG] csvPath (after trim) = "%s"', csvPath);

csvPath2 = fullfile(string(P.RESOURCE), char(csvPath));
logmsg(logf, '[MONTAGE][DEBUG] csvPath2 (P.RESOURCE + csv) = "%s"', string(csvPath2));

if isfolder(P.RESOURCE)
    D = dir(P.RESOURCE);
    names = string({D.name});
    logmsg(logf, '[MONTAGE][DEBUG] Files in P.RESOURCE = "%s"', strjoin(names, ", "));
else
    logmsg(logf, '[MONTAGE][DEBUG] P.RESOURCE folder does not exist.');
end

if ~isfile(csvPath)
    error('Montage CSV not found: %s', char(csvPath));
end

labs0 = string({EEG.chanlocs.labels});
logmsg(logf, '[DEBUG] nbchan before AB-select = %d', EEG.nbchan);
logmsg(logf, '[DEBUG] First 20 labels: %s', strjoin(labs0(1:min(20,end)), ", "));

wantA = "A" + string(1:32);
wantB = "B" + string(1:32);
present = upper(strtrim(labs0));

logmsg(logf, '[DEBUG] Present A chans: %s', strjoin(intersect(wantA, present), ", "));
logmsg(logf, '[DEBUG] Present B chans: %s', strjoin(intersect(wantB, present), ", "));

% Optionally select only A/B scalp channels
selectAB = true;
if isfield(cfg.exp.montage, 'select_ab_only')
    selectAB = logical(cfg.exp.montage.select_ab_only);
end

if selectAB
    labs = string({EEG.chanlocs.labels});
    labsU = upper(strtrim(labs));

    isAB = ~cellfun(@isempty, regexp(cellstr(labsU), '^(A|B)\d+$|^(A|B)0+\d+$', 'once'));
    keep = cellstr(labs(isAB));

    EEG = pop_select(EEG, 'channel', keep);
    EEG = eeg_checkset(EEG);
    logmsg(logf, '[MONTAGE] Selected A1-32 & B1-32 only. nbchan = %d', EEG.nbchan);
end

% Read mapping
T = readtable(char(csvPath), 'Delimiter', ',', 'TextType', 'string');
T.Properties.VariableNames = lower(string(T.Properties.VariableNames));
reqsCols = ["raw_label", "std_label"];
if ~all(ismember(reqsCols, string(T.Properties.VariableNames)))
    error('Montage CSV must contain columns: raw_label, std_label');
end

raw = string(T.raw_label);
std = string(T.std_label);

% Uniqueness checks
if numel(unique(raw)) ~= numel(raw)
    error('Montage CSV has duplicate raw_label entries.');
end
if numel(unique(std)) ~= numel(std)
    error('Montage CSV has duplicate std_label entires.');
end

% Apply relabel
cur = string({EEG.chanlocs.labels});
curU = upper(cur);
rawU = upper(raw);

missing = setdiff(rawU, curU);
if ~isempty(missing)
    error('Montage CSV expects channel not present in EEG: %s', strjoin(missing, ', '));
end

for i = 1:numel(rawU)
    idx = find(curU == rawU(i), 1);
    EEG.chanlocs(idx).labels = char(std(i));
end
EEG = eeg_checkset(EEG);

% Duplicate check (case-insensitive)
labs = string({EEG.chanlocs.labels});
if numel(unique(upper(labs))) ~= numel(labs)
    error('Relabel produced duplicate labels (case-insensitive). Check montage CSV.');
end

logmsg(logf, '[MONTAGE] Relabel complete from CSV: %s', char(csvPath));

% Optional ELP lookup
doLookup = true;
if isfield(cfg.exp.montage, 'do_lookup')
    doLookup = logical(cfg.exp.montage.do_lookup);
end

if doLookup
    elpPath = "";
    if isfield(P, 'CORE') && isfield(P.CORE, 'ELP_FILE')
        elpPath = string(P.CORE.ELP_FILE);
    end
    if strlength(elpPath) > 0 && exist(elpPath, 'file')
        EEG = pop_chanedit(EEG, 'lookup', char(elpPath));
        EEG = eeg_checkset(EEG);
        logmsg(logf, '[MONTAGE] Applied coord lookup from ELP: %s', elpPath);
    else
        logmsg(logf, '[WARN] ELP missing; skipping lookup (P.CORE.ELP_FILE).');
    end
end

write_channelmap_tsv(LOGS, subjid, EEG);
end