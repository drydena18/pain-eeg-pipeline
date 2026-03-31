function sessPaths = resolve_raw_session_files(P, cfg, subjid, nSessions, logf)
% RESOLVE_RAW_SESSION_FILES  Locate raw EEG files for every session
% V 1.0.0
%
% Tries three strategies in order for each session index (1..nSessions):
%
%   1. BIDS multi-session layout:
%        <INPUT.EXP>/sub-NNN/ses-NN/eeg/sub-NNN_ses-NN_task-TASK_eeg.bdf|eeg
%
%   2. cfg.preproc.concat.session_pattern — a printf format string whose
%      arguments are (subjid, session_index), e.g.:
%        "sub%03d/ses%02d/sub%03d_ses%02d_eeg.bdf"   -> (subjid, sess, subjid, sess)
%        "sub%03d_ses%02d.bdf"                        -> (subjid, sess)
%      The pattern is interpreted relative to P.INPUT.EXP.
%      The function tries formatting with 2 or 4 arguments automatically.
%
%   3. Recursive dir() fallback — finds all raw files for this subject
%      anywhere under P.INPUT.EXP, sorted alphabetically.  The k-th result
%      is used for session k.  This is the safest fallback for non-BIDS
%      datasets but requires that directory listing order matches session
%      order (usually true when filenames contain a session number).
%
% Inputs:
%   P          : paths struct from config_paths.m
%   cfg        : pipeline cfg struct
%   subjid     : integer subject ID
%   nSessions  : number of sessions to resolve (= cfg.preproc.concat.n_sessions)
%   logf       : (optional) log file handle
%
% Output:
%   sessPaths  : [nSessions x 1] cell array of char file paths.
%                Throws an error if any session file cannot be resolved.

if nargin < 5, logf = 1; end

sessPaths = cell(nSessions, 1);

subDir = sprintf('sub-%03d', subjid);

% Task name for BIDS candidate construction
task = "";
if isfield(cfg, 'exp') && isfield(cfg.exp, 'task') && strlength(string(cfg.exp.task)) > 0
    task = string(cfg.exp.task);
elseif isfield(P, 'EXP') && isfield(P.EXP, 'raw_dirname')
    task = string(P.EXP.raw_dirname);
else
    task = "task";
end

% Pattern override (optional)
pat = "";
if isfield(cfg, 'preproc') && isfield(cfg.preproc, 'concat') && ...
        isfield(cfg.preproc.concat, 'session_pattern')
    pat = string(cfg.preproc.concat.session_pattern);
end

% Raw extensions to try (ordered by preference)
exts = {'.bdf', '.BDF', '.eeg', '.EEG'};

% ----------------------------------------------------------------
% Fallback 3 pre-computation: collect all raw files for this subject
% ----------------------------------------------------------------
allRaw = {};
for e = 1:numel(exts)
    d = dir(fullfile(string(P.INPUT.EXP), subDir, '**', ['*' exts{e}]));
    if isempty(d)
        % Broaden to whole dataset
        d = dir(fullfile(string(P.INPUT.EXP), '**', [subDir '*' exts{e}]));
    end
    for k = 1:numel(d)
        allRaw{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
    end
end
allRaw = unique(allRaw);
allRaw = sort(allRaw);   % alphabetical -> session order when names contain sess index

% ----------------------------------------------------------------
% Resolve each session
% ----------------------------------------------------------------
for s = 1:nSessions
    sessStr2  = sprintf('ses-%02d', s);   % BIDS zero-padded (ses-01)
    sessStr1  = sprintf('ses-%d',   s);   % non-padded (ses-1)

    found = '';

    % ---- Strategy 1: BIDS layout ----
    for ext = exts
        bidsDir  = fullfile(string(P.INPUT.EXP), subDir, sessStr2, 'eeg');
        cands = {
            fullfile(bidsDir, sprintf('%s_%s_task-%s_eeg%s', subDir, sessStr2, task, ext{1}))
            fullfile(bidsDir, sprintf('%s_%s_eeg%s',          subDir, sessStr2,        ext{1}))
            % non-padded session variant
            fullfile(string(P.INPUT.EXP), subDir, sessStr1, 'eeg', ...
                sprintf('%s_%s_task-%s_eeg%s', subDir, sessStr1, task, ext{1}))
        };
        for c = 1:numel(cands)
            if exist(cands{c}, 'file')
                found = cands{c};
                break;
            end
        end
        if ~isempty(found), break; end
    end

    % ---- Strategy 2: cfg pattern ----
    if isempty(found) && strlength(pat) > 0
        try
            % Try with 4 args first (subjid, sess, subjid, sess)
            rel4 = sprintf(char(pat), subjid, s, subjid, s);
            cand4 = fullfile(string(P.INPUT.EXP), rel4);
            if exist(cand4, 'file')
                found = cand4;
            end
        catch
        end
        if isempty(found)
            try
                % Try with 2 args (subjid, sess)
                rel2 = sprintf(char(pat), subjid, s);
                cand2 = fullfile(string(P.INPUT.EXP), rel2);
                if exist(cand2, 'file')
                    found = cand2;
                end
            catch
            end
        end
    end

    % ---- Strategy 3: k-th file in sorted recursive listing ----
    if isempty(found) && s <= numel(allRaw)
        found = allRaw{s};
        logmsg(logf, '[CONCAT][WARN] sub-%03d sess %d: using fallback recursive file: %s', ...
            subjid, s, found);
    end

    if isempty(found)
        error('resolve_raw_session_files:NotFound', ...
            'sub-%03d session %d: no raw file found. Add cfg.preproc.concat.session_pattern or use BIDS layout.', ...
            subjid, s);
    end

    logmsg(logf, '[CONCAT] sub-%03d sess %d -> %s', subjid, s, found);
    sessPaths{s} = char(found);
end
end