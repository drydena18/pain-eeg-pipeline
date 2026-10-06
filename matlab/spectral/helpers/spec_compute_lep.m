function lep = spec_compute_lep(EEG, timesSec, F, lepCfg, outDir, subjid, logf)
% SPEC_COMPUTE_LEP  Save LEP waveforms + per trial x channel peak table
% V 2.0.0
%
% V2.0.0 changes vs V1.0.0:
%   - Peaks come from spec_lep_features with the SHARED windows
%     (cfg.features.lep: N2 150-350 ms, P2 250-500 ms), identical to the
%     source pipeline; column names now match the source CSV
%     (n2_amp, n2_lat_ms, p2_amp, p2_lat_ms, n2p2_amp, n2_mean, p2_mean).
%     The same columns are also written into sub-XXX_spectral_chan_by_trial.csv.
%   - cfg.spectral.lep.window_sec now only sets the stored waveform span.
%
% Saves in outDir:
%   sub-XXX_lep_trials.mat  data [nChan x nT x nTr], t_ms, chan_labels, fs
%   sub-XXX_lep_ga.mat      data [nChan x nT],       t_ms, chan_labels, fs
%   sub-XXX_lep_peaks.csv   one row per (trial, channel)

if nargin < 7, logf = 1; end

fs = EEG.srate;
mW = spec_window_mask(timesSec, lepCfg.window_sec(1), lepCfg.window_sec(2), fs);
chanLabels = spec_get_chanlabels(EEG);
nChan = EEG.nbchan;
nTr = EEG.trials;

trialData = double(EEG.data(:, mW, :));
gaData = mean(trialData, 3);
t_ms = timesSec(mW) * 1000;

lepTrial = struct('data', trialData, 't_ms', t_ms, 'chan_labels', {chanLabels}, 'fs', fs, 'subjid', subjid);
lepGA    = struct('data', gaData,    't_ms', t_ms, 'chan_labels', {chanLabels}, 'fs', fs, 'subjid', subjid);
save(fullfile(outDir, sprintf('sub-%03d_lep_trials.mat', subjid)), '-struct', 'lepTrial');
save(fullfile(outDir, sprintf('sub-%03d_lep_ga.mat', subjid)),     '-struct', 'lepGA');

lep = spec_lep_features(EEG.data, timesSec, fs, F.lep.n2_window_sec, F.lep.p2_window_sec);

nRows = nChan * nTr;
T = table();
T.subjid     = repmat(int32(subjid), nRows, 1);
T.trial      = repelem((1:nTr)', nChan);
T.chan_idx   = repmat((1:nChan)', nTr, 1);
T.chan_label = repmat(string(chanLabels(:)), nTr, 1);
fns = fieldnames(lep);
for k = 1:numel(fns)
    T.(fns{k}) = reshape(lep.(fns{k}), nRows, 1);
end
writetable(T, fullfile(outDir, sprintf('sub-%03d_lep_peaks.csv', subjid)));
spec_logmsg(logf, '[LEP] N2 %s s, P2 %s s; saved waveforms + peaks to %s', ...
    mat2str(F.lep.n2_window_sec(:)'), mat2str(F.lep.p2_window_sec(:)'), outDir);
end
