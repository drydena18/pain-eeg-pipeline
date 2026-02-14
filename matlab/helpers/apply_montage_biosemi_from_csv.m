function EEG = apply_montage_biosemi_from_csv(P, cfg, EEG, logf, subjid, LOGS)
% Expects csv columns: raw_label, std_label
% Writes audit file: LOGS/sub-###_channelmap_applied.tsv

% Resolve csv path
csvPath = "";
if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'csv')
    csvPath = string(cfg.exp.montage.csv);
end
if strlength(csvPath) == 0
    error('Montage enabled but cfg.exp.montage.csv is missing.');
end

% Prefer resource folder copy if relative
csvPath2 = fullfile(string(P.RESOURCE), char(csvPath));
if isfile(csvPath2)
    csvPath = string(csvPath2);
end
if ~isfile(csvPath)
    error('Montage CSV not found: %s', char(csvPath));
end

% Optionally select only A/B scalp channels
selectAB = true;
if isfield(cfg.exp.montage, 'select_ab_only')
    selectAB = logical(cfg.exp.montage.select_ab_only);
end

if selectAB
    keep = [
        arrayfun(@(x) sprintf('A%d', x), 1:32, 'UniformOutput', false), 
        arrayfun(@(x) sprintf('B%d', x), 1:32, 'UniformOutput', fale)];
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