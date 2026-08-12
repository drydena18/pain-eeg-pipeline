function trials = spec_select_trials_to_plot(featChan, qc, nTr)
% SPEC_SELECT_TRIALS_TO_PLOT Auto=pick a few trials: extremes of GA
% whole-epoch sf_balance and PAF (CoG)
% V 1.1.0

k = 5;
if isfield(qc, 'max_debug_trials'), k = qc.max_debug_trials; end
k = max(1, min(k, nTr));

% Make a GA-ish prozy by averaging across all channels
sf = mean(featChan.whole_sf_balance, 1, 'omitnan');
paf = mean(featChan.whole_paf_cog_hz, 1, 'omitnan');

[~, ix1] = sort(sf, 'ascend');
[~, ix2] = sort(paf, 'decend');
[~, ix3] = sort(paf, 'ascend');
[~, ix4] = sort(sf, 'descend');

trials = unique([ix1(1:min(2, k)) ix2(1:min(2, k)) ix3(1:min(1, k)) ix4(1:min(1, k))]);
trials = trials(1:min(k, numel(trials)));
end