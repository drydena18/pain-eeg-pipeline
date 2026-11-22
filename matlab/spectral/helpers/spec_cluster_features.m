function [featC, summC, auxC] = spec_cluster_features(X, timesSec, fs, F, aux, C)
% SPEC_CLUSTER_FEATURES  Features for electrode clusters (pooled channels)
% V 1.0.0
%
% Each cluster is treated as ONE unit with the same feature set and the
% same column names as a single channel, so its rows can be appended to the
% trial x channel and subject x channel CSVs (label = cluster name).
%
% Pooling rule (decided per feature type, never by averaging metric values):
%
%   Power features - pool POWER, then compute the metric.
%       pow_slow/fast/alpha : mean across the cluster's channels of each
%                             channel's filter-Hilbert window power
%       paf_cog_hz          : CoG of the cluster-mean Welch PSD
%       ratios / balance / psi_cog / ERD / p5_flag
%                           : from the pooled powers (spec_metrics_from_powers,
%                             spec_noise_stats and spec_erd_metrics on the
%                             cluster's own powers)
%     This is the rule already used for the whole-head GA unit in
%     spec_compute_subject_features. Averaging per-channel ratios instead
%     would give a different (and biased) number, because the ratio of means
%     is not the mean of ratios.
%
%   Waveform features - pool the SIGNAL (a "virtual electrode"), then compute.
%       slow-alpha phase, ITC, LEP (N2/P2)
%     Phase and evoked peaks are only defined for one time series. Neighbouring
%     electrodes over one hemisphere are close to in phase in the alpha band,
%     so the cluster-mean signal is a fair summary. Because power features
%     are pooled from each channel's own power, signal cancellation between
%     channels can only affect these waveform features, never the alpha
%     metrics that are screened.
%
%   Subject level (one row per cluster, as in the channel summary):
%       trial-mean pooled powers, CoG of the trial-mean pooled PSD, ratios from
%       those, ERD without p5, phase/LEP of the trial-mean virtual electrode,
%       ITC of the virtual electrode, TVI_alpha of the pooled pre_sf_balance
%       sequence. FOOOF is fitted by spectral_core on auxC.wholePsdTrialMean.
%
% Inputs:
%   X        : [nChan x nTime x nTr] epochs (EEG.data)
%   timesSec : [1 x nTime] s
%   fs       : Hz
%   F        : cfg.spectral.feat
%   aux      : 4th output of spec_compute_subject_features (needs f, pw, PS)
%   C        : output of spec_resolve_clusters
%
% Outputs:
%   featC : fields [nClus x nTr], same names as featChan
%   summC : fields [nClus x 1],  same names as summ (before FOOOF)
%   auxC  : names {1 x nClus}, ok [nClus x 1],
%           wholePsdTrialMean [nClus x nF] (NaN rows for clusters not ok)

if ~isfield(aux, 'pw') || ~isfield(aux, 'PS')
    error('spec_cluster_features:OldAux', ...
        'aux lacks pw/PS: spec_compute_subject_features must be V1.2.0 or later.');
end

timesSec = double(timesSec(:))';
nC  = numel(C);
nTr = size(X, 3);
nF  = numel(aux.f);

trialRows = cell(nC, 1);
subjRows  = cell(nC, 1);
psdRows   = nan(nC, nF);
ok        = false(nC, 1);

for c = 1:nC
    if ~C(c).ok
        continue;
    end
    [trialRows{c}, subjRows{c}, psdRows(c, :)] = local_one_cluster( ...
        double(X(C(c).idx, :, :)), C(c).idx, aux, timesSec, fs, F, nTr);
    ok(c) = true;
end

featC = local_stack(trialRows, ok, nC);
summC = local_stack(subjRows,  ok, nC);

auxC = struct();
auxC.names = {C.name};
auxC.ok = ok;
auxC.wholePsdTrialMean = psdRows;
end

% =========================================================================
function [trial, subj, psdRow] = local_one_cluster(Xs, idx, aux, timesSec, fs, F, nTr)
tb = F.band_power.trans_bw_hz;
wn = {'whole', 'pre', 'post'};

cpw = struct(); cP = struct(); cft = struct();
spw = struct(); sP = struct(); sft = struct();
for i = 1:numel(wn)
    w = wn{i};
    % Trial level: pool power across the cluster's channels -> [1 x nTr]
    cpw.(w) = struct( ...
        'slow',  mean(aux.pw.(w).slow(idx, :),  1, 'omitnan'), ...
        'fast',  mean(aux.pw.(w).fast(idx, :),  1, 'omitnan'), ...
        'alpha', mean(aux.pw.(w).alpha(idx, :), 1, 'omitnan'));
    cP.(w) = mean(aux.PS.(w)(idx, :, :), 1, 'omitnan');            % [1 x nF x nTr]
    cft.(w) = spec_metrics_from_powers(cpw.(w).slow, cpw.(w).fast, cpw.(w).alpha, ...
        spec_cog_from_psd(aux.f, cP.(w), F.bands.alpha_hz));

    % Subject level: trial mean of the pooled powers / PSD -> scalar, [1 x nF]
    spw.(w) = struct( ...
        'slow',  mean(cpw.(w).slow,  2, 'omitnan'), ...
        'fast',  mean(cpw.(w).fast,  2, 'omitnan'), ...
        'alpha', mean(cpw.(w).alpha, 2, 'omitnan'));
    sP.(w) = reshape(mean(cP.(w), 3, 'omitnan'), 1, []);
    sft.(w) = spec_metrics_from_powers(spw.(w).slow, spw.(w).fast, spw.(w).alpha, ...
        spec_cog_from_psd(aux.f, sP.(w), F.bands.alpha_hz));
end

noise = spec_noise_stats(cpw.pre.slow, cpw.pre.fast, F.erd.eps0_frac, F.erd.p5_percentile);

% Virtual electrode for waveform features
Xc = mean(Xs, 1, 'omitnan');                                         % [1 x nTime x nTr]
[ph, itc] = spec_slow_phase_features(Xc, timesSec, fs, F.bands.slow_hz, tb, ...
    F.phase.post_ref_sec, F.windows.post_sec);
lep = spec_lep_features(Xc, timesSec, fs, F.lep.n2_window_sec, F.lep.p2_window_sec);

trial = spec_merge_structs(struct(), ...
    spec_add_prefix(cft.whole, 'whole_'), ...
    spec_add_prefix(cft.pre,   'pre_'), ...
    spec_add_prefix(cft.post,  'post_'), ...
    spec_compute_metric_deltas(cft.pre, cft.post), ...
    spec_erd_metrics(cpw.pre, cpw.post, noise, true), ...
    ph, lep);

Xcm = mean(Xc, 3);
phS  = spec_slow_phase_features(Xcm, timesSec, fs, F.bands.slow_hz, tb, ...
    F.phase.post_ref_sec, F.windows.post_sec);
lepS = spec_lep_features(Xcm, timesSec, fs, F.lep.n2_window_sec, F.lep.p2_window_sec);
lepS.n_trials = nTr;
tvi = spec_compute_tvi_alpha(cft.pre.sf_balance, -1);

subj = spec_merge_structs(struct(), ...
    spec_add_prefix(sft.whole, 'whole_'), ...
    spec_add_prefix(sft.pre,   'pre_'), ...
    spec_add_prefix(sft.post,  'post_'), ...
    spec_compute_metric_deltas(sft.pre, sft.post), ...
    spec_erd_metrics(spw.pre, spw.post, noise, false), ...
    phS, lepS, itc, struct('TVI_alpha', tvi.tvi_alpha));

psdRow = sP.whole;
end

% =========================================================================
function S = local_stack(rows, ok, nC)
% Stack per-cluster structs of [1 x n] fields into [nClus x n]; clusters
% that were not computed get NaN rows.
S = struct();
first = find(ok, 1);
if isempty(first)
    return;
end
tmpl = rows{first};
fn = fieldnames(tmpl);
for k = 1:numel(fn)
    v = tmpl.(fn{k});
    M = nan(nC, numel(v));
    for c = find(ok(:))'
        M(c, :) = double(rows{c}.(fn{k})(:)');
    end
    S.(fn{k}) = M;
end
end