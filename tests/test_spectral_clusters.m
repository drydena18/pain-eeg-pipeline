function test_spectral_clusters()
% TEST_SPECTRAL_CLUSTERS  Synthetic-data checks for electrode clusters
% V 1.0.0  (MATLAB or Octave; no EEGLAB needed)
%
% Checks spec_resolve_clusters, spec_cluster_features (with
% spec_compute_subject_features V1.2.0) and spec_append_units:
%   - contra/ipsi naming follows the stimulated hand
%   - a 1-electrode cluster reproduces that channel exactly (every column)
%   - power pooled across electrodes; ratios from pooled powers
%   - known 9 Hz (left) / 11 Hz (right) alpha recovered; ERD and LEP sane
%   - 32-channel cap coverage, missing-electrode handling, name clashes
%
% Run from the repo root:  addpath('tests'); test_spectral_clusters

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'matlab', 'spectral', 'helpers'));
rng(1);
pass = 0; nChecks = 0;

labs64 = {'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz', ...
          'Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','CZ','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'};
fs = 500; t = -2.2:1/fs:2.0; nT = numel(t); nTr = 25; nCh = numel(labs64);
X = 0.5*randn(nCh, nT, nTr);
iL = find(ismember(upper(labs64), {'C3','CP1','CP5'}));
iR = find(ismember(upper(labs64), {'C4','CP2','CP6'}));
for tr = 1:nTr
  a = 1 + 0.3*randn; ph = 2*pi*rand;
  sL = a*sin(2*pi*9*t + ph);   sL(t > 0.2 & t < 0.8) = 0.4*sL(t > 0.2 & t < 0.8);  % slow alpha + ERD
  sR = a*sin(2*pi*11*t + ph);
  X(iL, :, tr) = X(iL, :, tr) + [1; 0.8; 0.6] .* sL;   % unequal weights -> ratio-of-means != mean-of-ratios
  X(iR, :, tr) = X(iR, :, tr) + [1; 0.8; 0.6] .* sR;
  lep = -3*exp(-((t-0.25)/0.03).^2) + 4*exp(-((t-0.4)/0.05).^2);
  X(:, :, tr) = X(:, :, tr) + lep;
end

F = struct();
F.windows    = struct('whole_sec', [-1.0, 2.0], 'pre_sec', [-1.0, -0.1], 'post_sec', [0.1, 0.8]);
F.bands      = struct('alpha_hz', [8, 12], 'slow_hz', [8, 10], 'fast_hz', [10, 12]);
F.band_power = struct('method', 'filter_hilbert', 'trans_bw_hz', 1.5);
F.psd        = struct('fmin_hz', 1, 'fmax_hz', 40, 'window_sec', 2.0, 'overlap_frac', 0.5, 'df_target_hz', 0.25);
F.erd        = struct('eps0_frac', 1e-3, 'p5_percentile', 5);
F.lep        = struct('n2_window_sec', [0.15, 0.35], 'p2_window_sec', [0.25, 0.50]);
F.phase      = struct('post_ref_sec', 0.2);

[featChan, featGA, summ, aux] = spec_compute_subject_features(X, t, fs, F);

% cfg as jsondecode would return it (column cellstrs), plus a 1-channel test cluster
cc = struct('enabled', true, 'stim_side', 'right', 'min_channels', 1);
cc.lateral = struct('S1', struct('left', {{'C3';'CP1';'CP5'}}, 'right', {{'C4';'CP2';'CP6'}}));
cc.fixed = struct('SM_furman', {{'C3';'Cz';'C4'}}, 'ONE_C3', {{'c3'}}, 'x_note', {{'ignored'}});
C = spec_resolve_clusters(cc, labs64);
names = {C.name};
[pass, nChecks] = tally(pass, nChecks, isequal(names, {'S1_contra','S1_ipsi','SM_furman','ONE_C3'}), 'unit names / order / comment key skipped');
[pass, nChecks] = tally(pass, nChecks, isequal(sort(C(1).idx(:))', sort(iL(:))'), 'stim right -> S1_contra = left hemisphere (C3/CP1/CP5)');
[pass, nChecks] = tally(pass, nChecks, isequal(sort(C(2).idx(:))', sort(iR(:))'), 'S1_ipsi = right hemisphere (C4/CP2/CP6)');
[pass, nChecks] = tally(pass, nChecks, any(strcmp(C(3).found, 'CZ')), 'case-insensitive match (Cz -> CZ)');

cc2 = cc; cc2.stim_side = 'left'; C2 = spec_resolve_clusters(cc2, labs64);
[pass, nChecks] = tally(pass, nChecks, isequal(sort(C2(1).idx(:))', sort(iR(:))'), 'stim left -> S1_contra = right hemisphere');

[featC, summC, auxC] = spec_cluster_features(X, t, fs, F, aux, C);
[pass, nChecks] = tally(pass, nChecks, isequal(sort(fieldnames(featC)), sort(fieldnames(featChan))), 'trial-level columns identical to channel columns');
[pass, nChecks] = tally(pass, nChecks, isequal(sort(fieldnames(summC)), sort(fieldnames(summ))), 'subject-level columns identical to channel summary');

% 1-electrode cluster must reproduce that channel exactly
c3 = find(strcmp(labs64, 'C3')); k = 4;
fn = fieldnames(featChan); mx = 0;
for i = 1:numel(fn), d = abs(featC.(fn{i})(k,:) - featChan.(fn{i})(c3,:)); mx = max([mx, d(~isnan(d))]); end
[pass, nChecks] = tally(pass, nChecks, mx < 1e-9, sprintf('1-electrode cluster == channel, all %d trial fields (max diff %.2e)', numel(fn), mx));
fn = fieldnames(summ); mx = 0;
for i = 1:numel(fn), d = abs(summC.(fn{i})(k) - summ.(fn{i})(c3)); mx = max([mx, d(~isnan(d))]); end
[pass, nChecks] = tally(pass, nChecks, mx < 1e-9, sprintf('1-electrode cluster == channel, all %d subject fields (max diff %.2e)', numel(fn), mx));

% Pooling rules
pm = mean(featChan.pre_pow_slow_alpha(iL, :), 1);
[pass, nChecks] = tally(pass, nChecks, max(abs(featC.pre_pow_slow_alpha(1,:) - pm)) < 1e-9, 'power = mean of channel powers');
ps = mean(featChan.pre_pow_slow_alpha(iL,:),1); pf = mean(featChan.pre_pow_fast_alpha(iL,:),1);
bal = (ps - pf) ./ (ps + pf + 1e-12);
[pass, nChecks] = tally(pass, nChecks, max(abs(featC.pre_sf_balance(1,:) - bal)) < 1e-9, 'sf_balance from pooled powers (ratio of means)');
[pass, nChecks] = tally(pass, nChecks, max(abs(featC.pre_sf_balance(1,:) - mean(featChan.pre_sf_balance(iL,:),1))) > 1e-6, '... which differs from mean of per-channel balances');
[pass, nChecks] = tally(pass, nChecks, mean(featC.pre_sf_balance(1,:)) > 0.3 && mean(featC.pre_sf_balance(2,:)) < -0.3, 'contra (9 Hz) slow-dominant, ipsi (11 Hz) fast-dominant');
[pass, nChecks] = tally(pass, nChecks, abs(mean(featC.pre_paf_cog_hz(1,:)) - 9) < 0.6 && abs(mean(featC.pre_paf_cog_hz(2,:)) - 11) < 0.6, ...
   sprintf('CoG contra %.2f Hz, ipsi %.2f Hz', mean(featC.pre_paf_cog_hz(1,:)), mean(featC.pre_paf_cog_hz(2,:))));
[pass, nChecks] = tally(pass, nChecks, mean(featC.erd_slow(1,:)) < -0.3, sprintf('contra slow-alpha ERD negative (%.2f)', mean(featC.erd_slow(1,:))));
[pass, nChecks] = tally(pass, nChecks, abs(mean(featC.n2_lat_ms(1,:)) - 250) < 15 && abs(mean(featC.p2_lat_ms(1,:)) - 400) < 20, 'LEP N2/P2 latencies from virtual electrode');
[pass, nChecks] = tally(pass, nChecks, all(size(auxC.wholePsdTrialMean) == [4, numel(aux.f)]), 'FOOOF input rows = clusters');

% spectral.units = "clusters": only member electrodes are processed
[members, Csub] = spec_cluster_subset(C, nCh);
expMembers = unique([iL(:); iR(:); find(strcmp(labs64, 'CZ'))]);
[pass, nChecks] = tally(pass, nChecks, isequal(members(:), expMembers(:)), sprintf('clusters-only processes %d of %d electrodes', numel(members), nCh));
tFull = tic; [~, ~, ~, auxF] = spec_compute_subject_features(X, t, fs, F); tFull = toc(tFull);
tSub  = tic; [~, ~, ~, auxS] = spec_compute_subject_features(X(members, :, :), t, fs, F);
[fcS, scS] = spec_cluster_features(X(members, :, :), t, fs, F, auxS, Csub); tSub = toc(tSub);
fn = fieldnames(featC); mx = 0;
for i = 1:numel(fn), d = abs(fcS.(fn{i}) - featC.(fn{i})); mx = max([mx; d(~isnan(d))]); end
[pass, nChecks] = tally(pass, nChecks, mx < 1e-10, sprintf('clusters-only == all-electrode run, trial level (max diff %.1e)', mx));
fn = fieldnames(summC); mx = 0;
for i = 1:numel(fn), d = abs(scS.(fn{i}) - summC.(fn{i})); mx = max([mx; d(~isnan(d))]); end
[pass, nChecks] = tally(pass, nChecks, mx < 1e-10, sprintf('clusters-only == all-electrode run, subject level (max diff %.1e)', mx));
fprintf('      feature extraction: all %d electrodes %.2f s | clusters only %.2f s (%.1fx faster)\n', nCh, tFull, tSub, tFull / tSub);

% Append to channel tables
[S, L] = spec_append_units(featChan, labs64, featC, names);
[pass, nChecks] = tally(pass, nChecks, numel(L) == nCh + 4 && strcmp(L{end-3}, 'S1_contra'), 'labels appended');
[pass, nChecks] = tally(pass, nChecks, isequal(S.pre_sf_balance(1:nCh,:), featChan.pre_sf_balance) && isequal(S.pre_sf_balance(nCh+1:end,:), featC.pre_sf_balance), 'channel rows unchanged, cluster rows appended');
featC2 = rmfield(featC, 'n2_amp'); S2 = spec_append_units(featChan, labs64, featC2, names);
[pass, nChecks] = tally(pass, nChecks, all(all(isnan(S2.n2_amp(nCh+1:end,:)))), 'field missing for clusters -> NaN rows');
[Ss] = spec_append_units(summ, labs64, summC, names);
[pass, nChecks] = tally(pass, nChecks, size(Ss.TVI_alpha,1) == nCh + 4, 'subject table appended');

% 32-channel cap: CP3 requested, not present; min_channels enforcement
labs32 = {'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','POz','O1','Oz','O2'};
cc3 = cc; cc3.min_channels = 2; cc3.fixed = struct('DENSE', {{'C1';'CP3';'C3'}});
C3 = spec_resolve_clusters(cc3, labs32);
[pass, nChecks] = tally(pass, nChecks, all([C3(1:2).ok]) && numel(C3(1).idx) == 3, 'S1 clusters complete on a 32-channel 10-20 cap');
[pass, nChecks] = tally(pass, nChecks, ~C3(3).ok && isequal(C3(3).missing, {'C1','CP3'}), 'missing electrodes reported; below min_channels -> not ok');
X32 = 0.5*randn(32, nT, 6);
[~, ~, ~, a32] = spec_compute_subject_features(X32, t, fs, F);
[fc32, sc32, ax32] = spec_cluster_features(X32, t, fs, F, a32, C3);
[pass, nChecks] = tally(pass, nChecks, all(isnan(fc32.pre_sf_balance(3,:))) && all(~isnan(fc32.pre_sf_balance(1,:))), 'not-ok cluster -> NaN rows, others computed');

% Name clash
cc4 = cc; cc4.fixed = struct('C3', {{'C3';'C1'}});
try, spec_resolve_clusters(cc4, labs64); r = 0; catch, r = 1; end
[pass, nChecks] = tally(pass, nChecks, r == 1, 'cluster named like a real channel -> error');
fprintf('\n%d / %d checks passed\n', pass, nChecks);
if pass ~= nChecks
    error('test_spectral_clusters:Failed', '%d check(s) failed.', nChecks - pass);
end
end

function [pass, n] = tally(pass, n, cond, msg)
n = n + 1;
if cond
    fprintf('PASS  %s\n', msg);
    pass = pass + 1;
else
    fprintf('FAIL  %s\n', msg);
end
end