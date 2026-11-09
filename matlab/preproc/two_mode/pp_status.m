function [state, info] = pp_status(cfg, subj)
%PP_STATUS Furthest VALID state for one participant.
%
%   raw           no committed ICA-ready file
%   corrupt       marker exists but the files don't match it (info.reason says why)
%   ica_ready     ICA-ready file valid; review packet missing or stale
%   review_ready  ICLabel + figures ready; waiting for manual rejection
%   reviewed      IC decisions saved for THIS decomposition; Mode 2b pending
%   final         final dataset committed for current decisions + cfg.mode2b
%
%   info.ica_stale is true when the ICA-ready file is intact but was made
%   with a different cfg.mode1. It is still used downstream (so reviewed
%   ICs aren't thrown away) unless pp_run is called with 'rerun_stale'.

P = pp_paths(cfg, subj);
info = struct('ica', [], 'ica_stale', false, 'review', [], 'decisions', [], 'final', [], 'reason', '');

[icaState, info.ica, info.reason] = local_ica(P, pp_hash(cfg.mode1));
if any(strcmp(icaState, {'raw', 'corrupt'}))
    state = icaState;
    return;
end
info.ica_stale = strcmp(icaState, 'stale');
w = info.ica.weights_hash;

% Decisions are bound to the decomposition (weights fingerprint), not to Mode 2a.
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
function [st, m, why] = local_ica(P, h1)
st = 'raw'; m = []; why = '';
if ~isfile(P.ica_marker), return; end

st = 'corrupt';
setf = fullfile(P.ica_dir, P.ica_set);
fdtf = [setf(1:end-4) '.fdt'];
if ~isfile(setf), why = 'ICA-ready .set missing'; return; end
if ~isfile(fdtf), why = 'ICA-ready .fdt missing'; return; end

try
    mk = pp_read_json(P.ica_marker);
catch ME
    why = ['marker unreadable: ' ME.message]; return;
end
try
    H = local_header(setf);       % reads the .set directly: no eeg_checkset side effects
catch ME
    why = ['.set header unreadable: ' ME.message]; return;
end

d = dir(fdtf);
expected = 4 * double(H.nbchan) * double(H.pnts) * double(H.trials);   % float32
if d.bytes ~= expected
    why = sprintf('.fdt is %d bytes, expected %d (truncated?)', d.bytes, expected); return;
end
if isempty(H.icaweights), why = 'no ICA weights in .set'; return; end
if size(H.icaweights, 1) ~= mk.rank
    why = sprintf('%d components in .set, marker says rank %d', size(H.icaweights, 1), mk.rank); return;
end
if ~isequal(double(H.icachansind(:)'), double(mk.icachansind(:)'))
    why = 'icachansind differs from marker'; return;
end

fp = pp_ica_fingerprint(H.icaweights, H.icasphere);
if isfield(mk, 'fp_version') && mk.fp_version >= 2
    if ~strcmp(fp, mk.weights_hash)
        why = 'ICA weights in .set do not match the marker'; return;
    end
else
    % Marker from the first version, which hashed the raw weights before
    % pop_saveset rescaled them, so it can never match. The structural
    % checks above passed: adopt the scale-invariant fingerprint.
    mk.weights_hash_legacy = mk.weights_hash;
    mk.weights_hash = fp;
    mk.fp_version = 2;
    pp_write_json(P.ica_marker, mk);
    fprintf('  %s: ICA-ready marker upgraded to scale-invariant fingerprint.\n', mk.subject);
end

m = mk;
if strcmp(mk.cfg_hash, h1), st = 'ok'; else, st = 'stale'; end
end

function H = local_header(setf)
need = {'nbchan', 'pnts', 'trials', 'icaweights', 'icasphere', 'icachansind'};
vars = who('-file', setf);
if ismember('EEG', vars)              % saved as a single EEG variable
    S = load(setf, '-mat', 'EEG');
    H = S.EEG;
else                                  % saved with -struct (EEGLAB default)
    S = load(setf, '-mat', need{ismember(need, vars)});
    H = S;
end
miss = need(~isfield(H, need));
if ~isempty(miss), error('missing fields: %s', strjoin(miss, ', ')); end
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
