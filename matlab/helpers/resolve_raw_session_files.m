function [sessPaths, sessSidecars] = resolve_raw_session_files(P, cfg, subjid, nSessions, logf)
% RESOLVE_RAW_SESSION_FILES  Locate raw EEG files and BIDS sidecars for every session
% V 1.1.0
%
% Tries three strategies in order for each session index (1..nSessions):
%
%   1. BIDS multi-session layout:
%        <INPUT.EXP>/sub-NNN/ses-NN/eeg/sub-NNN_ses-NN_task-TASK_eeg.set|bdf|eeg
%
%   2. cfg.preproc.concat.session_pattern — a printf format string whose
%      arguments are (subjid, session_index), e.g.:
%        "sub%03d/ses%02d/sub%03d_ses%02d_eeg.set"   -> (subjid, sess, subjid, sess)
%        "sub%03d_ses%02d.set"                        -> (subjid, sess)
%      The pattern is interpreted relative to P.INPUT.EXP.
%      The function tries formatting with 2 or 4 arguments automatically.
%
%   3. Recursive dir() fallback — finds all raw files for this subject
%      anywhere under P.INPUT.EXP, sorted alphabetically.  The k-th result
%      is used for session k.  This is the safest fallback for non-BIDS
%      datasets but requires that directory listing order matches session
%      order (usually true when filenames contain a session number).
%
% After resolving each raw file, the function also searches for BIDS sidecar
% files in the same eeg/ directory:
%   channels.tsv   — channel labels (applied to EEG.chanlocs before .elp lookup)
%   coordsystem.json — logged for provenance only; not parsed
%
% Inputs:
%   P          : paths struct from config_paths.m
%   cfg        : pipeline cfg struct
%   subjid     : integer subject ID
%   nSessions  : number of sessions to resolve (= cfg.preproc.concat.n_sessions)
%   logf       : (optional) log file handle
%
% Outputs:
%   sessPaths    : [nSessions x 1] cell array of char raw file paths.
%                  Throws an error if any session file cannot be resolved.
%   sessSidecars : [nSessions x 1] cell array of structs with fields:
%                    .channels_tsv     — char path, or '' if not found
%                    .coordsystem_json — char path, or '' if not found

if nargin < 5, logf = 1; end

sessPaths    = cell(nSessions, 1);
sessSidecars = cell(nSessions, 1);

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
% .set first: EEGLAB native format (also handles .fdt automatically)
exts = {'.set', '.bdf', '.BDF', '.eeg', '.EEG'};

% ----------------------------------------------------------------
% Fallback 3 pre-computation: collect all raw files for this subject.
% Break on the first extension that yields results so we never mix
% .set and .eeg files from the same sessions in allRaw.
% ----------------------------------------------------------------
allRaw = {};
for e = 1:numel(exts)
    d = dir(fullfile(string(P.INPUT.EXP), subDir, '**', ['*' exts{e}]));
    if isempty(d)
        % Broaden to whole dataset
        d = dir(fullfile(string(P.INPUT.EXP), '**', [subDir '*' exts{e}]));
    end
    if ~isempty(d)
        for k = 1:numel(d)
            allRaw{end+1} = fullfile(d(k).folder, d(k).name); %#ok<AGROW>
        end
        break;  % stop at first extension that finds files
    end
end
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
            % non-padded session variant (ses-1 instead of ses-01)
            fullfile(string(P.INPUT.EXP), subDir, sessStr1, 'eeg', ...
                sprintf('%s_%s_task-%s_eeg%s', subDir, sessStr1, task, ext{1}))
            fullfile(string(P.INPUT.EXP), subDir, sessStr1, 'eeg', ...
                sprintf('%s_%s_eeg%s', subDir, sessStr1, ext{1}))
        };
        for c = 1:numel(cands)
            if exist(char(cands{c}), 'file')
                found = char(cands{c});
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
            if exist(char(cand4), 'file')
                found = char(cand4);
            end
        catch
        end
        if isempty(found)
            try
                % Try with 2 args (subjid, sess)
                rel2 = sprintf(char(pat), subjid, s);
                cand2 = fullfile(string(P.INPUT.EXP), rel2);
                if exist(char(cand2), 'file')
                    found = char(cand2);
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

    % ---- Sidecar resolution ----
    % Derive the eeg/ directory from the resolved raw file path.
    % Both channels.tsv and coordsystem.json follow the same BIDS naming
    % convention rooted in the same folder.
    eegDir   = fileparts(found);
    % Strip the raw file suffix to get the BIDS basename stem
    [~, rawName] = fileparts(found);   % e.g. sub-001_ses-1_task-29ByANT_eeg
    % Replace trailing _eeg (from _eeg.set) to build sidecar stems
    baseStem = regexprep(rawName, '_eeg$', '');  % sub-001_ses-1_task-29ByANT

    chanTsv  = fullfile(eegDir, [baseStem '_channels.tsv']);
    coordJson = fullfile(eegDir, [baseStem '_coordsystem.json']);

    sidecar = struct();
    if exist(char(chanTsv), 'file')
        sidecar.channels_tsv = char(chanTsv);
        logmsg(logf, '[CONCAT] sub-%03d sess %d: channels.tsv -> %s', subjid, s, chanTsv);
    else
        sidecar.channels_tsv = '';
        logmsg(logf, '[CONCAT][WARN] sub-%03d sess %d: channels.tsv not found: %s', subjid, s, chanTsv);
    end

    if exist(char(coordJson), 'file')
        sidecar.coordsystem_json = char(coordJson);
        logmsg(logf, '[CONCAT] sub-%03d sess %d: coordsystem.json -> %s', subjid, s, coordJson);
    else
        sidecar.coordsystem_json = '';
    end

    sessSidecars{s} = sidecar;
end
end