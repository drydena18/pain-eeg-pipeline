function [featChan, featGA, summ, aux] = spec_compute_subject_features(X, timesSec, fs, F)
% SPEC_COMPUTE_SUBJECT_FEATURES  All channel-space features for one subject
% V 1.1.0  (V1.1.0: whole window from F.windows.whole_sec, not the full epoch)
%
% Python twin: python/source/source_core.compute_subject_features. Kept free
% of I/O so tests/parity_test.py can compare both pipelines on identical data.
% Every definition below is the same as in source space, with unit =
% channel instead of ROI.
%
% Inputs:
%   X        : [nChan x nTime x nTr] baseline-corrected epochs (EEG.data)
%   timesSec : [1 x nTime] seconds (EEG.times / 1000)
%   fs       : Hz
%   F        : cfg.spectral.feat (shared features)
%
% Outputs (unprefixed families get whole_/pre_/post_ prefixes):
%   featChan : trial x channel, fields [nChan x nTr]
%              whole_/pre_/post_<10 families>, delta_<10>, erd_slow, erd_fast,
%              erd_pow_alpha_total, erd_asym, p5_flag, slow_phase, sin_phase,
%              cos_phase, slow_phase_post, sin_phase_post, cos_phase_post,
%              n2_amp, n2_lat_ms, p2_amp, p2_lat_ms, n2p2_amp, n2_mean, p2_mean
%              -> identical names to sub-XXX_source_trial.csv
%   featGA   : trial x spatial grand average (channel mean), fields [1 x nTr]
%              powers = channel mean of channel powers, CoG from the
%              channel-mean PSD, ratios from the mean powers; ERD/p5 on the
%              GA unit's own powers; phase = circular mean across channels.
%              (No LEP: the channel mean of average-referenced data is 0.)
%   summ     : subject x channel, fields [nChan x 1]
%              powers = trial mean, CoG from the trial-mean PSD, ratios from
%              mean powers, ERD (no p5), phase/LEP of the trial-mean waveform,
%              ITC, TVI_alpha -> identical names to sub-XXX_source_ga.csv
%              (FOOOF added by spectral_core)
%   aux      : f, wholePsdTrialMean [nChan x nF], wholePsd [nChan x nF x nTr],
%              phaseT0 [nChan x nTr] (for r_cl), gaWholePsd [nF x nTr]

timesSec = double(timesSec(:))';
X = double(X);
nCh = size(X, 1);
nTr = size(X, 3);
tb = F.band_power.trans_bw_hz;

bands = struct('slow', F.bands.slow_hz, 'fast', F.bands.fast_hz, 'alpha', F.bands.alpha_hz);
PT = spec_band_power_tcs(X, fs, bands, tb);

wins = struct('whole', F.windows.whole_sec, 'pre', F.windows.pre_sec, 'post', F.windows.post_sec);
wn = fieldnames(wins);
ft = struct(); pw = struct(); PS = struct();
for i = 1:numel(wn)
    [ft.(wn{i}), pw.(wn{i}), f, PS.(wn{i})] = spec_window_alpha_features(PT, X, timesSec, wins.(wn{i}), fs, F);
end

% ------------------------------------------------------------------
% Trial x channel
% ------------------------------------------------------------------
noise = spec_noise_stats(pw.pre.slow, pw.pre.fast, F.erd.eps0_frac, F.erd.p5_percentile);
[ph, itc] = spec_slow_phase_features(X, timesSec, fs, F.bands.slow_hz, tb, F.phase.post_ref_sec, F.windows.post_sec);
lep = spec_lep_features(X, timesSec, fs, F.lep.n2_window_sec, F.lep.p2_window_sec);

featChan = spec_merge_structs(struct(), ...
    spec_add_prefix(ft.whole, 'whole_'), ...
    spec_add_prefix(ft.pre,   'pre_'), ...
    spec_add_prefix(ft.post,  'post_'), ...
    spec_compute_metric_deltas(ft.pre, ft.post), ...
    spec_erd_metrics(pw.pre, pw.post, noise, true), ...
    ph, lep);

% ------------------------------------------------------------------
% Trial x spatial GA (channel mean as one unit)
% ------------------------------------------------------------------
gft = struct(); gpw = struct();
for i = 1:numel(wn)
    w = wn{i};
    gpw.(w) = struct( ...
        'slow',  mean(pw.(w).slow,  1, 'omitnan'), ...
        'fast',  mean(pw.(w).fast,  1, 'omitnan'), ...
        'alpha', mean(pw.(w).alpha, 1, 'omitnan'));
    gP = mean(PS.(w), 1, 'omitnan');                       % [1 x nF x nTr]
    gft.(w) = spec_metrics_from_powers(gpw.(w).slow, gpw.(w).fast, gpw.(w).alpha, ...
        spec_cog_from_psd(f, gP, F.bands.alpha_hz));
end
gNoise = spec_noise_stats(gpw.pre.slow, gpw.pre.fast, F.erd.eps0_frac, F.erd.p5_percentile);
gPh = struct();
phn = fieldnames(ph);
for i = 1:numel(phn)
    if startsWith(phn{i}, 'slow_phase')
        cm = angle(mean(exp(1i * ph.(phn{i})), 1, 'omitnan'));
        suffix = phn{i}(numel('slow_phase') + 1:end);
        gPh.(['slow_phase' suffix]) = cm;
        gPh.(['sin_phase' suffix])  = sin(cm);
        gPh.(['cos_phase' suffix])  = cos(cm);
    end
end
featGA = spec_merge_structs(struct(), ...
    spec_add_prefix(gft.whole, 'whole_'), ...
    spec_add_prefix(gft.pre,   'pre_'), ...
    spec_add_prefix(gft.post,  'post_'), ...
    spec_compute_metric_deltas(gft.pre, gft.post), ...
    spec_erd_metrics(gpw.pre, gpw.post, gNoise, true), ...
    gPh);

% ------------------------------------------------------------------
% Subject x channel (trial mean) - parity with source GA (subject x ROI)
% ------------------------------------------------------------------
sft = struct(); spw = struct(); sP = struct();
for i = 1:numel(wn)
    w = wn{i};
    spw.(w) = struct( ...
        'slow',  mean(pw.(w).slow,  2, 'omitnan'), ...
        'fast',  mean(pw.(w).fast,  2, 'omitnan'), ...
        'alpha', mean(pw.(w).alpha, 2, 'omitnan'));
    sP.(w) = mean(PS.(w), 3, 'omitnan');                   % [nCh x nF]
    sft.(w) = spec_metrics_from_powers(spw.(w).slow, spw.(w).fast, spw.(w).alpha, ...
        spec_cog_from_psd(f, sP.(w), F.bands.alpha_hz));
end
Xm = mean(X, 3);
phS = spec_slow_phase_features(Xm, timesSec, fs, F.bands.slow_hz, tb, F.phase.post_ref_sec, F.windows.post_sec);
lepS = spec_lep_features(Xm, timesSec, fs, F.lep.n2_window_sec, F.lep.p2_window_sec);
lepS.n_trials = repmat(nTr, nCh, 1);

tviVal = nan(nCh, 1);
for c = 1:nCh
    t = spec_compute_tvi_alpha(ft.pre.sf_balance(c, :), -1);
    tviVal(c) = t.tvi_alpha;
end

summ = spec_merge_structs(struct(), ...
    spec_add_prefix(sft.whole, 'whole_'), ...
    spec_add_prefix(sft.pre,   'pre_'), ...
    spec_add_prefix(sft.post,  'post_'), ...
    spec_compute_metric_deltas(sft.pre, sft.post), ...
    spec_erd_metrics(spw.pre, spw.post, noise, false), ...
    phS, lepS, itc, struct('TVI_alpha', tviVal));

aux = struct();
aux.f = f;
aux.wholePsd = PS.whole;
aux.wholePsdTrialMean = sP.whole;
aux.gaWholePsd = reshape(mean(PS.whole, 1, 'omitnan'), numel(f), nTr);
aux.phaseT0 = ph.slow_phase;
aux.gaPreSfBalance = gft.pre.sf_balance(:);
end
