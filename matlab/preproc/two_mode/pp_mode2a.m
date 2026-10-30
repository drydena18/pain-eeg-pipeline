function meta = pp_mode2a(cfg, subj)
%PP_MODE2A Automated review prep: ICLabel, suggested rejections, QC figures.
%   Writes to the review folder (and any QC folders the prep hook uses);
%   the ICA-ready file is never modified.

P = pp_paths(cfg, subj);
c = cfg.mode2a;
pp_ensure_dirs(cfg, subj);
[~, info] = pp_status(cfg, subj);

rtc = [];
if isfield(cfg, 'runtime_fn') && ~isempty(cfg.runtime_fn)
    [cfg.rt, rtc] = cfg.runtime_fn(cfg, subj, 'mode2a');
end %#ok<NASGU>

EEG = pop_loadset('filename', P.ica_set, 'filepath', P.ica_dir);

pqc = struct();
if isfield(c, 'prep_fn') && ~isempty(c.prep_fn)
    [EEG, suggested, pqc] = c.prep_fn(EEG, cfg, subj);       % e.g. st_icqc (CNED helpers)
else
    EEG = pop_iclabel(EEG, c.iclabel_version);
    L0  = EEG.etc.ic_classification.ICLabel;
    suggested = local_flag(L0.classifications, L0.classes, c.flag);
end
suggested = double(suggested(:)');

L = EEG.etc.ic_classification.ICLabel;
probs   = L.classifications;
classes = L.classes;

nShow = min(c.n_show, size(EEG.icaweights, 1));
pp_plot_ic_grid(EEG, nShow, probs, classes, suggested, P.ic_png, subj);
reasons = {};
if isfield(pqc, 'reasons'), reasons = pqc.reasons; end
save(P.iclabel_mat, 'probs', 'classes', 'suggested', 'reasons', '-v7');

meta = struct();
meta.subject      = subj;
meta.stage        = 'review_ready';
meta.cfg_hash     = pp_hash(c);
meta.weights_hash = info.ica.weights_hash;
meta.n_ic         = size(probs, 1);
meta.n_show       = nShow;
meta.suggested    = suggested;
meta.reasons      = reasons;
meta.finished     = char(datetime('now'), 'yyyy-MM-dd''T''HH:mm:ss');
pp_write_json(P.review_marker, meta);   % commit
end

function idx = local_flag(probs, classes, flag)
names = regexprep(classes, '\s', '');
hit = false(size(probs, 1), 1);
f = fieldnames(flag);
for i = 1:numel(f)
    j = find(strcmpi(names, f{i}), 1);
    if isempty(j)
        error('pp:flag', 'Unknown ICLabel class "%s" (have: %s).', f{i}, strjoin(names, ', '));
    end
    hit = hit | probs(:, j) >= flag.(f{i});
end
idx = find(hit)';
end