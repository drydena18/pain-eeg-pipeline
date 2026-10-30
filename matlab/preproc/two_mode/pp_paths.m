function P = pp_paths(cfg, subj)
%PP_PATHS File locations for one participant, inside the per-subject folder
%   preproc_core uses:  <out>/sub-XXX/<folder>/<prefix><NNN>_...
if isempty(cfg.paths.out), error('pp:paths', 'Set cfg.paths.out.'); end
s = cfg.paths.sub;

pre = '';
if isfield(cfg.paths, 'prefix'), pre = char(cfg.paths.prefix); end
num = regexp(subj, '\d+', 'match', 'once');
if isempty(num)
    stem = [pre subj];
else
    stem = sprintf('%s%03d', pre, str2double(num));    % same as local_fname
end

subDir = fullfile(cfg.paths.out, subj);
P.subject_dir = subDir;
P.stem        = stem;
P.dirs.ica    = fullfile(subDir, s.ica);
P.dirs.review = fullfile(subDir, s.review);
P.dirs.final  = fullfile(subDir, s.final);
P.dirs.logs   = fullfile(cfg.paths.out, s.logs);       % experiment level

P.ica_dir       = P.dirs.ica;
P.ica_set       = [stem '_icaready.set'];               % name only (pop_saveset style)
P.ica_marker    = fullfile(P.dirs.ica, [stem '_icaready.json']);

P.iclabel_mat   = fullfile(P.dirs.review, [stem '_iclabel.mat']);
P.ic_png        = fullfile(P.dirs.review, [stem '_ic_topos.png']);
P.review_marker = fullfile(P.dirs.review, [stem '_review.json']);
P.decisions     = fullfile(P.dirs.review, [stem '_icdecisions.json']);

% Mode 2b's final file is whatever its last stage saved (08_base for the
% CNED pipeline); the marker records its path. final_set is the fallback.
P.final_dir     = P.dirs.final;
P.final_set     = [stem '_final.set'];
P.final_marker  = fullfile(P.dirs.final, [stem '_pp_final.json']);

P.log           = fullfile(P.dirs.logs, 'pp_errors.log');
P.manifest      = fullfile(P.dirs.logs, 'pp_manifest.csv');
end