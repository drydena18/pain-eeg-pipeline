function pp_review(cfg, varargin)
%PP_REVIEW Manual IC rejection for every participant in 'review_ready'.
%   Shows the topoplot grid made in Mode 2a, loads the ICA-ready dataset,
%   and calls cfg.review.ui_fn. Each decision is saved immediately, so you
%   can stop at any point and pick up later. When the session ends, Mode 2b
%   runs for everyone with saved decisions (cfg.review.autorun_2b).
%
%   ui_fn signature:
%     [rejected, ok] = fn(subj, probs, classes, suggested, i, n, EEG, cfg)
%   ok = false pauses the session (nothing saved for this participant).

p = inputParser;
p.addParameter('subjects', cfg.subjects);
p.addParameter('run_2b', cfg.review.autorun_2b);
p.parse(varargin{:});
o = p.Results;

subs = cellstr(o.subjects);
todo = subs(cellfun(@(s) strcmp(pp_status(cfg, s), 'review_ready'), subs));
fprintf('%d participant(s) awaiting IC review.\n', numel(todo));

for i = 1:numel(todo)
    subj = todo{i};
    P = pp_paths(cfg, subj);
    [~, info] = pp_status(cfg, subj);
    R = load(P.iclabel_mat);
    EEG = pop_loadset('filename', P.ica_set, 'filepath', P.ica_dir);

    fig = figure('Name', sprintf('%s - IC topographies', subj), 'NumberTitle', 'off');
    imshow(imread(P.ic_png), 'Border', 'tight');
    [rej, ok] = cfg.review.ui_fn(subj, R.probs, R.classes, R.suggested, i, numel(todo), EEG, cfg);
    if isgraphics(fig), close(fig); end

    if ~ok
        fprintf('Review paused at %s. Saved decisions are kept; run the review again to continue.\n', subj);
        break;
    end

    dec = struct();
    dec.subject      = subj;
    dec.weights_hash = info.ica.weights_hash;   % binds decisions to this decomposition
    dec.n_ic         = size(R.probs, 1);
    dec.rejected     = sort(double(rej(:)))';
    dec.suggested    = double(R.suggested(:))';
    dec.added        = setdiff(dec.rejected, dec.suggested);
    dec.overridden   = setdiff(dec.suggested, dec.rejected);
    dec.reviewer     = cfg.review.reviewer;
    dec.reviewed     = char(datetime('now'), 'yyyy-MM-dd''T''HH:mm:ss');
    pp_write_json(P.decisions, dec);
    fprintf('  %s: rejected %d IC(s).\n', subj, numel(dec.rejected));
end

if o.run_2b
    pp_run(cfg, 'modes', {'2b'}, 'subjects', subs);
end
end