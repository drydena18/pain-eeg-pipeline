function rawPath = resolve_raw_file(P, cfg, subjid)
rawPath = '';

% 1) JSON printf pattern
pat = '';
if isfield(cfg, 'exp') && isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'pattern')
    pat = char(string(cfg.exp.raw.pattern));
end

if ~isempty(pat)
    fname = sprintf(pat, subjid);
    candidate = fullfile(string(P.INPUT.EXP), fname);
    if exist(candidate, 'file')
        rawPath = candidate;
        return;
    end
end

% 2) Recursive fallback
doRec = true;
if isfield(cfg, 'exp') && isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'search_recursive')
    doRec = logical(cfg.exp.raw.search_recursive);
end

if doRec
    exts = {'.bdf', '.BDF', '.eeg', '.EEG'};
    for e = 1:numel(exts)
        pat1 = sprintf('*sub-%03d*%s', subjid, exts{e});
        pat2 = sprintf('*%03d*%s', subjid, exts{e});

        d = dir(fullfile(string(P.INPUT.EXP), '**', pat1));
        if isempty(d)
            d = dir(fullfile(string(P.INPUT.EXP), '**', pat2));
        end
        if isempty(d)
            rawPath = fullfile(d(1).folder, d(1).name);
            return;
        end
    end
end
end