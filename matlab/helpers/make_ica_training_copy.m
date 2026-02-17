function [EEGtrain, segInfo] = make_ica_training_copy(EEG, cfg, logf)
EEGtrain = EEG;
segInfo = struct('rempoved', false, 'n_intervals', 0, 'pct_time', 0, 'intervals', []);

if ~isfield(cfg, 'preproc') || ~isfield(cfg.preproc, 'initrej') || ~isfield(cfg.preproc.initrej, 'badseg') || ~cfg.preproc.initrej.badseg.enabled
    logmsg(logf, '[ICA-TRAIN] badseg disabled; using full data.');
    return;
end

thr = cfg.preproc.initrej.badseg.threshold_uv;
logmsg(logf, '[ICA-TRAIN] Detecting bad segments (thr = %.1f uV) for training copy.', thr);

x = double(EEG.data);
badSamp = any(abs(x) > thr, 1);
intervals = mask_to_intervals(badSamp)';

segInfo.intervals = intervals;
segInfo.n_intervals = size(intervals, 1);
segInfo.pct_time = 100 * (sum(badSamp) / numel(badSamp));

if isempty(intervals)
    logmsg(logf, '[ICA-TRAIN] No bad segments detected; using full data.');
    return;
end

logmsg(logf, '[ICA-TRAIN] Detected %d intervals (%.2f%% of samples).', segInfo.n_intervals, segInfo.pct_time);

doRemove = prompt_yesno('Remove detected bad segments from ICA training copy? (y/n) [n]:', false);
if ~doRemove
    logmsg(logf, '[ICA-TRAIN] Keeping all segments (manual decision).');
    return;
end

EEGtrain = pop_select(EEGtrain, 'nopoint', intervals);
EEGtrain = eeg_checkset(EEGtrain);

segInfo.removed = true;
logmsg(logf, '[ICA-TRAIN] Removed bad segments from training copy.');
end