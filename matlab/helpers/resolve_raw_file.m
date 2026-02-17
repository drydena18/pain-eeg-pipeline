function rawPath = resolve_raw_file(P, cfg, subjid)
rawPath = "";

subDir = sprintf('sub-%03d', subjid);

% ----------------
% 0) BIDS first
% ----------------
% Prefer cfg.exp.task if available; else infer from registry raw_dirname
task = "";
if isfield(cfg, 'exp') && isfield(cfg.exp, 'task') && strlength(string(cfg.exp.task)) > 0
    task = string(cfg.exp.task);
elseif isfield(P, 'EXP') && isfield(P.EXP, 'raw_dirname')
    task = string(P.EXP.raw_dirname);
else
    task = "task";
end

bidsEEGDir = fullfile(string(P.INPUT.EXP), subDir, "eeg");

cand = [
    fullfile(bidsEEGDir, sprintf('%s_task-%s_eeg.bdf', subDir, task))
    fullfile(bidsEEGDir, sprintf('%s_task-%s_eeg.BDF', subDir, task))
    fullfile(bidsEEGDir, sprintf('%s_task-%s_eeg.eeg', subDir, task))
    fullfile(bidsEEGDir, sprintf('%s_task-%s_eeg.EEG'))
    ];

for i = 1:numel(cand)
    if exist(cand{i}, 'file')
        rawPath = string(cand{i});
        return;
    end
end

% --------------------------------
% 1) Optional JSON printf pattern
% --------------------------------
pat = "";
if isfield(cfg, 'exp') && isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'pattern')
    pat = string(cfg.exp.raw.pattern);
end

% If the pattern is relative, interpret relative to P.INPUT.EXP
if strlength(pat) > 0
    % Allow patterns that include directories
    try
        rel = sprintf(char(pat), subjid);
        cand2 = fullfile(string(P.INPUT.EXP), rel);
        if exist(cand2, 'file')
            rawPath = string(cand2);
            return;
        end
    catch
        % ignore pattern errors, fall through
    end
end

% -----------------------------
% 2) Recursive fallback search
% -----------------------------
doRec = true;
if isfield(cfg, 'exp') && isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'search_recursive')
    doRec = logical(cfg.exp.raw.search_recursive);
end
if ~doRec
    return;
end

exts = {'.bdf', '.BDF', '.eeg', '.EEG'};
for e = 1:numel(exts)
    % Search within expected subject folder first (faster)
    d = dir(fullfile(string(P.INPUT.EXP), subDir, '**', ['*' exts{e}]));
    if ~isempty(d)
        rawPath = string(fullfile(d(1).folder, d(1).name));
        return;
    end
end

% If still not found, broaden search across entire dataset (slow)
for e = 1:numel(exts)
    pat1 = sprintf('%s*%s', subDir, exts{e});
    d = dir(fullfile(string(P.INPUT.EXP), '**', pat1));
    if ~isempty(d)
        rawPath = string(fullfile(d(1).folder, d(1).name));
        return;
    end
end

rawPath = char(rawPath);

end