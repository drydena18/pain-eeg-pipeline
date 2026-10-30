function [removedICs, ok] = st_review_ui(subj, ~, ~, suggested, i, n, EEG, cfg)
%ST_REVIEW_UI The manual half of preproc_core's ICLabel block, run at review
%   time: IC activation scroll (ica.scroll) -> prompt_ic_reject (Enter =
%   remove none) -> post-removal channel scroll (ica.confirm_scroll).
%   A y/n prompt first lets you stop the session; earlier decisions are kept.
ica = cfg.exp.preproc.ica;
id  = str2double(regexp(subj, '\d+', 'match', 'once'));
nIC = size(EEG.icaweights, 1);

fprintf('\n===== IC review %s (%d/%d) | %d ICs | suggested: %s =====\n', subj, i, n, nIC, vec2str(suggested));
fprintf('QC packets: %s\n', fullfile(cfg.P.SUBJECT_DIR(id), 'QC', sprintf('sub-%03d_icqc', id)));

removedICs = [];
ok = prompt_yesno(sprintf('Review %s now? (y = review, n = stop for now) [y]: ', subj), true);
if ~ok, return; end

prompt_scroll_eeg(EEG, local_block(ica, 'scroll'), 1, 'components', 'ICA component scroll');

while true
    removedICs = prompt_ic_reject(suggested);
    bad = removedICs(removedICs < 1 | removedICs > nIC | removedICs ~= round(removedICs));
    if isempty(bad), break; end
    fprintf(2, 'Out of range (1-%d): %s. Try again.\n', nIC, vec2str(bad));
end

block = local_block(ica, 'confirm_scroll');
if isfield(block, 'enabled') && logical(block.enabled)
    tmp = EEG;
    if ~isempty(removedICs), tmp = pop_subcomp(EEG, removedICs, 0); end
    prompt_scroll_eeg(tmp, block, 1, 'channels', 'ICA post-reject confirmation scroll');
end
end

function b = local_block(ica, name)
b = struct('enabled', false);
if isfield(ica, name), b = ica.(name); end
end