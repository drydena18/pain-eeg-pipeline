function meta = pp_mode2b(cfg, subj)
%PP_MODE2B Apply saved IC decisions, run the post-ICA stages, commit.
%   Stage hooks save their own stage files (06_ica, 07_epoch, 08_base, ...);
%   the last file a stage reports in qc.saved is the final dataset. It is
%   re-saved with the run metadata in EEG.etc.pipeline_final, then the
%   marker is written. Never resumes from old stage files: outputs depend
%   on the IC decisions, which may have changed.

P = pp_paths(cfg, subj);
c = cfg.mode2b;
pp_ensure_dirs(cfg, subj);
[~, info] = pp_status(cfg, subj);
dec = info.decisions;
t0  = datetime('now');

pp_delete_if(P.final_marker);

rtc = [];
if isfield(cfg, 'runtime_fn') && ~isempty(cfg.runtime_fn)
    [cfg.rt, rtc] = cfg.runtime_fn(cfg, subj, 'mode2b');
    cfg.rt.resume = false;
end %#ok<NASGU>

EEG = pop_loadset('filename', P.ica_set, 'filepath', P.ica_dir);
if ~strcmp(pp_hash([EEG.icaweights(:); EEG.icasphere(:)]), dec.weights_hash)
    error('pp:decisions', 'IC decisions do not match the ICA weights on disk.');
end
if ~isfield(EEG.etc, 'pp_tags') && isfield(info.ica, 'tags')
    EEG.etc.pp_tags = local_cell(info.ica.tags);
end

rej = double(dec.rejected(:)');
nIC = size(EEG.icaweights, 1);
if any(rej < 1 | rej > nIC), error('pp:decisions', 'Rejected IC index out of range.'); end
if ~isempty(rej)
    EEG = pop_subcomp(EEG, rej, 0);
    EEG = eeg_checkset(EEG);
end
EEG.etc.pp_ic_rejected  = rej;
EEG.etc.pp_ic_suggested = double(dec.suggested(:)');

qc = struct();
finalPath = '';
for k = 1:numel(c.stages)
    st = c.stages(k);
    [EEG, sqc] = st.fn(EEG, cfg, subj);
    if isempty(sqc), sqc = struct(); end
    if isfield(sqc, 'saved') && ~isempty(sqc.saved), finalPath = char(sqc.saved); end
    qc.stages.(st.name) = sqc;
end
if isempty(finalPath)
    finalPath = fullfile(P.final_dir, P.final_set);
end
[fdir, fname, fext] = fileparts(finalPath);

meta = struct();
meta.subject        = subj;
meta.stage          = 'final';
meta.cfg_hash       = pp_hash(c);
meta.cfg_hash_mode1 = info.ica.cfg_hash;
meta.decisions_hash = pp_hash(fileread(P.decisions));
meta.mode2b_version = c.version;
meta.final_set_path = finalPath;
meta.n_ic           = nIC;
meta.rejected       = rej;
meta.n_trials       = EEG.trials;
meta.started        = char(t0, 'yyyy-MM-dd''T''HH:mm:ss');
meta.finished       = char(datetime('now'), 'yyyy-MM-dd''T''HH:mm:ss');
meta.env            = pp_env(cfg);
meta.qc             = qc;

EEG.etc.pipeline_final            = meta;   % Mode 1 metadata stays in EEG.etc.pipeline
EEG.etc.pipeline_final.decisions  = dec;
EEG.etc.pipeline_final.cfg_mode2b = pp_sanitize(c);
pop_saveset(EEG, 'filename', [fname fext], 'filepath', fdir, 'savemode', 'twofiles');

mj = meta;
try
    jsonencode(pp_sanitize(mj.qc));
catch
    mj.qc = 'not JSON-encodable; see EEG.etc.pipeline_final.qc';
end
pp_write_json(P.final_marker, mj);   % commit
end

function c = local_cell(x)
if isempty(x), c = {}; elseif ischar(x), c = {x}; else, c = cellstr(x(:)'); end
end