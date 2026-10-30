function meta = pp_mode1(cfg, subj)
%PP_MODE1 Raw -> ICA-ready dataset (weights transferred from clean copy).
%   The JSON (JavaScript Object Notation) marker is written LAST: it is the
%   commit. A crash anywhere before it leaves the participant in 'raw' or
%   'corrupt', and the next pp_run redoes it. Stage hooks may resume from
%   their own stage files (00-05) unless the ICA was stale (cfg.mode1 changed).

P  = pp_paths(cfg, subj);
c  = cfg.mode1;
t0 = datetime('now');
pp_ensure_dirs(cfg, subj);

[~, info0] = pp_status(cfg, subj);
resume = ~info0.ica_stale;     % config changed -> don't reuse old stage files

% Un-commit everything that depends on the old decomposition.
pp_delete_if(P.ica_marker);
pp_delete_if(P.review_marker);
pp_delete_if(P.final_marker);
pp_quarantine(P.decisions);    % kept on disk, renamed *_stale_<time>.json

rtc = [];
if isfield(cfg, 'runtime_fn') && ~isempty(cfg.runtime_fn)
    [cfg.rt, rtc] = cfg.runtime_fn(cfg, subj, 'mode1');
    cfg.rt.resume = resume;
end %#ok<NASGU> rtc closes the log on exit

% ---- Preprocessing stages ------------------------------------------------
EEG = c.load_fn(cfg, subj);
chanlocs_raw = EEG.chanlocs;
qc = struct('n_chan_raw', EEG.nbchan, 'srate_raw', EEG.srate, 'pnts_raw', EEG.pnts);
bad = {}; interp = {}; reported = false; interpUnknown = false;

for k = 1:numel(c.stages)
    st = c.stages(k);
    [EEG, sqc] = st.fn(EEG, cfg, subj);
    if isempty(sqc), sqc = struct(); end
    if isfield(sqc, 'bad_labels'),    bad    = [bad,    local_cell(sqc.bad_labels)];  end %#ok<AGROW>
    if isfield(sqc, 'interp_labels'), interp = [interp, local_cell(sqc.interp_labels)]; reported = true; end %#ok<AGROW>
    if isfield(sqc, 'interp_unknown') && sqc.interp_unknown, interpUnknown = true; end
    qc.stages.(st.name) = local_strip(sqc);
end
if ~reported && ~interpUnknown
    error('pp:rank', ['No stage returned qc.interp_labels. The bad-channel stage must ' ...
        'report it ({} if none) - the ICA rank depends on it.']);
end
bad = unique(bad, 'stable'); interp = unique(interp, 'stable');

% ---- Rank ------------------------------------------------------------------
% Reref BEFORE interpolation: the average-reference constraint includes the
% bad channels, so replacing them with interpolated signals breaks it and
% the rank loss is absorbed by the interpolated channels:
%   nInterp > 0  ->  rank = nChan - nInterp
%   nInterp = 0  ->  rank = nChan - 1
nInterp = numel(interp);
if interpUnknown
    rBook = NaN;   % reused stage file that didn't record its interpolation
else
    absorbed = c.avg_ref && c.reref_before_interp && nInterp > 0;
    rBook    = EEG.nbchan - nInterp - double(c.avg_ref && ~absorbed);
end
rEig = pp_eig_rank(EEG.data, c.rank_eig_tol);

% ---- Clean copy + rank-aware extended Infomax ------------------------------
EEGc = c.clean_copy_fn(EEG, cfg, subj);
if EEGc.nbchan ~= EEG.nbchan || ~isequal({EEGc.chanlocs.labels}, {EEG.chanlocs.labels})
    error('pp:cleancopy', 'Clean copy must keep the same channels in the same order.');
end
if isfield(EEGc, 'etc') && isstruct(EEGc.etc) && isfield(EEGc.etc, 'pp_train_seg')
    qc.clean_copy = rmfield(EEGc.etc.pp_train_seg, intersect({'intervals'}, fieldnames(EEGc.etc.pp_train_seg)));
end
rEigC = pp_eig_rank(EEGc.data, c.rank_eig_tol);
r = min([rBook, rEig, rEigC]);           % min ignores NaN
if isnan(rBook)
    warning('pp:rank', '%s: interpolation not recorded in reused 05_initrej file; using eigenvalue rank %d.', subj, r);
elseif rBook ~= rEig || rBook ~= rEigC
    warning('pp:rank', '%s: bookkeeping rank %d, eigenvalue rank %d (data) / %d (clean copy); using %d.', ...
        subj, rBook, rEig, rEigC, r);
end

seed = c.ica_seed;
if ischar(seed) || isstring(seed)          % 'subject' -> rng(subject number), as preproc_core
    seed = str2double(regexp(subj, '\d+', 'match', 'once'));
end
rng(seed, 'twister');
args = {'icatype', 'runica', 'extended', c.ica_extended};
if r < EEGc.nbchan, args = [args, {'pca', r}]; end
local_log(cfg, '[PP-ICA] rank %d (bookkeeping %g, eig %d / clean copy %d), seed %d', r, rBook, rEig, rEigC, seed);
EEGc = pop_runica(EEGc, args{:});

if ~isequal(EEGc.icachansind(:)', 1:EEG.nbchan)
    error('pp:cleancopy', 'Unexpected icachansind on clean copy; cannot transfer weights.');
end

% ---- Transfer weights to the main dataset ----------------------------------
EEG.icaweights  = EEGc.icaweights;
EEG.icasphere   = EEGc.icasphere;
EEG.icachansind = EEGc.icachansind;
EEG.icawinv = []; EEG.icaact = [];
savedChanlocs = EEG.chanlocs;
EEG = eeg_checkset(EEG, 'ica');            % recomputes icawinv from weights * sphere
if isempty(EEG.chanlocs) && ~isempty(savedChanlocs)
    EEG.chanlocs = savedChanlocs;          % same guard as preproc_core
end
EEG.icaact = [];
clear EEGc

% ---- Metadata --------------------------------------------------------------
tags = {};
if isfield(EEG, 'etc') && isstruct(EEG.etc) && isfield(EEG.etc, 'pp_tags'), tags = EEG.etc.pp_tags; end

meta = struct();
meta.subject            = subj;
meta.stage              = 'ica_ready';
meta.cfg_hash           = pp_hash(c);
meta.mode1_version      = c.version;
meta.rank               = size(EEG.icaweights, 1);
meta.rank_bookkeeping   = rBook;
meta.rank_eig           = rEig;
meta.rank_eig_cleancopy = rEigC;
meta.icachansind        = EEG.icachansind(:)';
meta.weights_hash       = pp_hash([EEG.icaweights(:); EEG.icasphere(:)]);
meta.n_chan             = EEG.nbchan;
meta.chan_labels        = {EEG.chanlocs.labels};
meta.bad_labels         = bad;
meta.interp_labels      = interp;
meta.interp_unknown     = interpUnknown;
meta.tags               = tags;
meta.ica_extended       = c.ica_extended;
meta.ica_seed           = seed;
meta.srate              = EEG.srate;
meta.pnts               = EEG.pnts;
meta.started            = char(t0, 'yyyy-MM-dd''T''HH:mm:ss');
meta.finished           = char(datetime('now'), 'yyyy-MM-dd''T''HH:mm:ss');
meta.env                = pp_env(cfg);
meta.qc                 = qc;

EEG.etc.pipeline              = meta;
EEG.etc.pipeline.chanlocs_raw = chanlocs_raw;
EEG.etc.pipeline.cfg_mode1    = pp_sanitize(c);

pop_saveset(EEG, 'filename', P.ica_set, 'filepath', P.ica_dir, 'savemode', 'twofiles');
local_log(cfg, '[PP-SAVE] ICA-ready: %s', fullfile(P.ica_dir, P.ica_set));

% Commit.
mj = meta;
try
    jsonencode(pp_sanitize(mj.qc));
catch
    mj.qc = 'not JSON-encodable; see EEG.etc.pipeline.qc';
end
pp_write_json(P.ica_marker, mj);
end

% -------------------------------------------------------------------------
function c = local_cell(x)
if isempty(x),      c = {};
elseif ischar(x),   c = {x};
elseif isstring(x), c = cellstr(x(:))';
else,               c = x(:)';
end
end

function s = local_strip(s)
% Keep QC small: drop large numeric arrays (they stay in EEG.etc where stages put them).
f = fieldnames(s);
for i = 1:numel(f)
    v = s.(f{i});
    if isnumeric(v) && numel(v) > 64, s = rmfield(s, f{i}); end
end
end

function local_log(cfg, fmt, varargin)
fprintf([fmt '\n'], varargin{:});
if isfield(cfg, 'rt') && isfield(cfg.rt, 'logf') && cfg.rt.logf > 2
    logmsg(cfg.rt.logf, fmt, varargin{:});
end
end