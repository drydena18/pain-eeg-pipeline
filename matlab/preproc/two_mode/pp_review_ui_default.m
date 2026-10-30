function [rejected, ok] = pp_review_ui_default(subj, probs, classes, suggested, i, n, varargin)
%PP_REVIEW_UI_DEFAULT Generic selection dialog (the CNED pipeline uses
%   st_review_ui instead, set in pp_hooks.m).
%   [rejected, ok] = fn(subj, probs, classes, suggested, i, n, EEG, cfg)

[p, j] = max(probs, [], 2);
items = arrayfun(@(k) sprintf('IC %3d   %-13s %.2f', k, classes{j(k)}, p(k)), ...
    (1:size(probs, 1))', 'UniformOutput', false);

args = {'ListString', items, 'SelectionMode', 'multiple', ...
        'Name', sprintf('%s  (%d/%d)', subj, i, n), ...
        'PromptString', {'Select ICs to REJECT (suggestions pre-ticked).', 'Cancel = pause review.'}, ...
        'ListSize', [340 520]};
if ~isempty(suggested), args = [args, {'InitialValue', suggested}]; end

[rejected, ok] = listdlg(args{:});
ok = logical(ok);
end