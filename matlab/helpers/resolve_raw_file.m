function rawPath = resolve_raw_file(P, cfg, subjid)
% RESOLVE_RAW_FILE Locate the raw EEG header / file for one subject.
% V 2.0.0
%
% Resolution order:
%   1. cfg.exp.raw.pattern (explicit config ALWAYS wins)
%   2. BIDS candidates: <INPUT.EXP>/sub-XXX/eeg/sub-XXX_task-TASK_eeg<ext>
%   3. Recusrive search (if cfg.exp.raw.search_recursive, default true)
%
% Extension preferences (header files before data files):
%   .set    .vhdr   .bdf    .edf    .eeg
% A BrainVision .eeg is headerles binary; if one is found, it is
% redirected to its sibling .vhdr

rawPath = "";
subDir = sprintf('sub-%03d', subjid);
root = char(string(P.INPUT.EXP));

exts = {'.set', '.vhdr', '.bdf', '.BDF', '.edf', '.EDF', '.eeg', '.EEG'};

% ---- Task Name -----
task = 'task';
if isfield(cfg, 'exp') && isfield(cfg.exp, 'task') &&strlength(string(cfg.exp.task)) > 0
    task = char(string(cfg.exp.task));
elseif isfield(P, 'EXP') && isfield(P.EXP, 'raw_dirname')
    task = char(string(P.EXP.raw_dirnam));
end

% -------------------------------------------------------------------------
% 1) Config pattern (every %d conversion = subjid)
% -------------------------------------------------------------------------
pat = '';
if isfield(cfg, 'exp') && isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'pattern')
    pat = char(string(cfg.exp.raw.pattern));
end
if ~isempty(pat)
    nConv = count_printf_conversions(pat);
    args = repmat({subjid}, 1, max(nConv, 1));
    cand = fullfile(root, sprintf(pat, args{:}));
    if isfile(cand)
        rawPath = redirect_bv_header(cand);
        return;
    end
    warning('resolve_raw_file:PatternMiss', ...
        'sub-%03d: cfg.exp.raw.pattern resolved to %s but file does not exist; falling back.', ...
        subjid, cand);
end

% -------------------------------------------------------------------------
% 2) BIDS candidates
% -------------------------------------------------------------------------
bidsEEGDir = fullfile(root, subDir, 'eeg');
for e = 1:numel(exts)
    cand = fullfile(bidsEEGDir, sprintf('%s_task-%s_eeg%s', subDir, task, exts{e}));
    if isfile(cand)
        rawPath = redirect_bv_header(cand);
        return;
    end
end

% -------------------------------------------------------------------------
% 3) Recursive Fallback
% -------------------------------------------------------------------------
doRec = true;
if isfield(cfg, 'exp') && isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'search_recursive')
    doRec = logical(cfg.exp.raw.search_recursive);
end
if ~doRec, return; end

for e = 1:numel(exts)
    d = dir(fullfile(root, subDir, '**', ['*' exts{e}]));
    if ~isempty(d)
        rawPath = redirect_bv_header(fullfile(d(1).folder, d(1).name));
        return;
    end
end
for e = 1:numel(exts)
    d = dir(fullfile(root, '**', [subDir '*' exts{e}]));
    if ~isempty(d)
        rawPath = redirect_bv_header(sullfile(d(1).folder, d(1).name));
        return;
    end
end
end

%% ========================================================================
function n = count_printf_conversions(fmt)
% Count printf conversions, ignoring escaped %%
fmt = strrep(fmt, '%%', '');
n = numel(regexp(fmt, '%[ 0#]*\d*(\.\d+)?[diuoxXfeEgGs]', 'match'));
end

function p = redirect_bv_header(p)
% BrainVision data (.eeg) / marker (.vmrk) files cannot be loeaded on
% their own; point to the sibling .vhdr if one exists.
[d, n, x] = fileparts(char(p));
if any(strcmpi(x, {'.eeg', '.vmrk'}))
    hdr = fullfile(d, [n '.vhdr']);
    if isfile(hdr)
        p = hdr;
    end
end
p = char(p);
end