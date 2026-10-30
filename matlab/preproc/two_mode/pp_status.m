function [state, info] = pp_status(cfg, subj)
%PP_STATUS Furthest VALID state for one participant.
%
%   raw       no committed ICA-ready file
%   corrupt   marker exists but the files don't match it (re-run Mode 1)
%   ica_ready ICA-ready file valid; review packet missing or stale
%   review_ready  ICLabel + figures ready; waiting for manual rejection
%   reviewed  IC decisions saved for THIS decomposition; Mode 2b pending
%   final     final dataset committed for current decisions + cfg.mode2b
%
%   info.ica_stale is true when the ICA-ready file is intact but was made
%   with a different cfg.mode1. It is still used downstream (so reviewed
%   ICs aren't thrown away) unless pp_run is called with 'rerun_stale'.

P = pp_paths(cfg, subj);
info = struct('ica', [], 'ica_stale', false, 'review', [], 'decisions', [], 'final', []);

[icaState, info.ica] = local_ica(P, pp_hash(cfg.mode1));
if any(strcmp(icaState, {'raw', 'corrupt'}))
    state = icaState;
    return;
end
info.ica_stale = strcmp(icaState, 'stale');
w = info.ica.weights_hash;

% Decisions are bound to the decomposition (weights hash), not to Mode 2a.
if isfile(P.decisions)
    d = pp_read_json(P.decisions);
    if strcmp(d.weights_hash, w)
        info.decisions = d;
        if local_final_ok(P, info.ica.cfg_hash, pp_hash(cfg.mode2b), pp_hash(fileread(P.decisions)))
            info.final = pp_read_json(P.final_marker);
            state = 'final';
        else
            state = 'reviewed';
        end
        return;
    end
end

if isfile(P.review_marker) && isfile(P.ic_png) && isfile(P.iclabel_mat)
    r = pp_read_json(P.review_marker);
    if strcmp(r.weights_hash, w) && strcmp(r.cfg_hash, pp_hash(cfg.mode2a))
        info.review = r;
        state = 'review_ready';
        return;
    end
end

state = 'ica_ready';
end

% -------------------------------------------------------------------------
function [st, m] = local_ica(P, h1)
st = 'raw'; m = [];
if ~isfile(P.ica_marker), return; end

st = 'corrupt';
setf = fullfile(P.ica_dir, P.ica_set);
fdtf = [setf(1:end-4) '.fdt'];
if ~(isfile(setf) && isfile(fdtf)), return; end

try
    mk = pp_read_json(P.ica_marker);
    [~, EEG] = evalc('pop_loadset(''filename'', P.ica_set, ''filepath'', P.ica_dir, ''loadmode'', ''info'')');
catch
    return;
end

d = dir(fdtf);
ok = d.bytes == 4 * EEG.nbchan * EEG.pnts * EEG.trials ...          % float32, not truncated
  && ~isempty(EEG.icaweights) ...
  && size(EEG.icaweights, 1) == mk.rank ...
  && isequal(double(EEG.icachansind(:)'), double(mk.icachansind(:)')) ...
  && strcmp(pp_hash([EEG.icaweights(:); EEG.icasphere(:)]), mk.weights_hash);
if ~ok, return; end

m = mk;
if strcmp(mk.cfg_hash, h1), st = 'ok'; else, st = 'stale'; end
end

function ok = local_final_ok(P, h1, h2b, hdec)
ok = false;
if ~isfile(P.final_marker), return; end
m = pp_read_json(P.final_marker);
setf = fullfile(P.final_dir, P.final_set);
if isfield(m, 'final_set_path') && ~isempty(m.final_set_path), setf = m.final_set_path; end
if ~(isfile(setf) && isfile([setf(1:end-4) '.fdt'])), return; end
ok = strcmp(m.cfg_hash_mode1, h1) && strcmp(m.cfg_hash, h2b) && strcmp(m.decisions_hash, hdec);
end