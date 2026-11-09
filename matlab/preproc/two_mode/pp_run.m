function T = pp_run(cfg, varargin)
%PP_RUN Unattended runner. Advances every participant as far as it can
%   without a human: Mode 1 -> Mode 2a, and Mode 2b for anyone whose IC
%   decisions already exist. Safe to re-run at any time; completed work
%   is detected from commit markers and skipped.
%
%   pp_run(cfg)                          % everything that can run
%   pp_run(cfg, 'modes', {'2b'})         % only apply finished reviews
%   pp_run(cfg, 'subjects', {'sub-017'})
%   pp_run(cfg, 'rerun_stale', true)     % redo ICA where cfg.mode1 changed
%
%   Returns the manifest table (also written to logs/pp_manifest.csv).

p = inputParser;
p.addParameter('modes', {'1', '2a', '2b'});
p.addParameter('subjects', cfg.subjects);
p.addParameter('rerun_stale', false);
p.parse(varargin{:});
o = p.Results;

subs = cellstr(o.subjects);
pp_ensure_dirs(cfg);
errs  = repmat({''}, numel(subs), 1);
order = {'1', '2a', '2b'};

for si = 1:numel(subs)
    subj = subs{si};
    for mi = 1:numel(order)
        mode = order{mi};
        if ~ismember(mode, o.modes), continue; end
        [state, info] = pp_status(cfg, subj);
        if ~local_needs(mode, state, info, o.rerun_stale), continue; end

        fprintf('[%s] Mode %s: %s\n', char(datetime('now', 'Format', 'HH:mm:ss')), mode, subj);
        try
            switch mode
                case '1',  pp_mode1(cfg, subj);
                case '2a', pp_mode2a(cfg, subj);
                case '2b', pp_mode2b(cfg, subj);
            end
        catch ME
            errs{si} = sprintf('Mode %s: %s%s', mode, ME.message, local_where(ME));
            pp_log_error(cfg, subj, mode, ME);
            fprintf(2, '  FAILED %s (Mode %s): %s\n', subj, mode, ME.message);
            break;   % don't run later modes for this participant
        end
    end
end

T = pp_manifest(cfg, subs, errs);
local_summary(T, cfg);
end

% -------------------------------------------------------------------------
function tf = local_needs(mode, state, info, rerunStale)
switch mode
    case '1',  tf = any(strcmp(state, {'raw', 'corrupt'})) || (rerunStale && info.ica_stale);
    case '2a', tf = strcmp(state, 'ica_ready');
    case '2b', tf = strcmp(state, 'reviewed');
    otherwise, error('pp:mode', 'Unknown mode "%s".', mode);
end
end

function local_summary(T, cfg)
fprintf('\n---- Summary ----\n');
states = {'final', 'reviewed', 'review_ready', 'ica_ready', 'raw', 'corrupt'};
for k = 1:numel(states)
    n = sum(strcmp(T.status, states{k}));
    if n, fprintf('  %-13s %d\n', states{k}, n); end
end
failed = T.subject(~cellfun(@isempty, T.last_error));
if ~isempty(failed)
    fprintf(2, '  failed: %s (see logs/pp_errors.log)\n', strjoin(failed', ', '));
end
stale = T.subject(T.ica_stale);
if ~isempty(stale)
    fprintf('  stale ICA (cfg.mode1 changed; kept, not re-run): %s\n', strjoin(stale', ', '));
end
nReview = sum(strcmp(T.status, 'review_ready'));
if nReview
    if isfield(cfg, 'exp_id')
        fprintf('\nNext: preproc_2mode(''%s'', ''review'')   %d participant(s) ready for IC rejection.\n', cfg.exp_id, nReview);
    else
        fprintf('\nNext: pp_review(cfg)   %d participant(s) ready for IC rejection.\n', nReview);
    end
end
end

function w = local_where(ME)
% First stack frame outside MATLAB/EEGLAB internals, as ' [file:line]'.
w = '';
for k = 1:numel(ME.stack)
    f = ME.stack(k).file;
    if contains(f, matlabroot) || contains(lower(f), 'eeglab'), continue; end
    [~, n, e] = fileparts(f);
    w = sprintf(' [%s%s:%d]', n, e, ME.stack(k).line);
    return;
end
if ~isempty(ME.stack)
    [~, n, e] = fileparts(ME.stack(1).file);
    w = sprintf(' [%s%s:%d]', n, e, ME.stack(1).line);
end
end
