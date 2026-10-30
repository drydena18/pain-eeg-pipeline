function EEGtrain = st_clean_copy(EEG, cfg, ~)
%ST_CLEAN_COPY make_ica_training_copy without the prompts. Same detection
%   (|x| > preproc.initrej.badseg.threshold_uv on any channel, intervals via
%   mask_to_intervals, removed with pop_select 'nopoint'); the train scroll
%   and prompt_yesno are replaced by cfg.mode1.badseg_policy:
%     'remove'  drop the detected segments from the training copy
%     'keep'    train on all data (the prompt's default)
%   Segment info is returned in EEGtrain.etc.pp_train_seg.
pre    = cfg.exp.preproc;
logf   = cfg.rt.logf;
policy = lower(char(cfg.mode1.badseg_policy));

EEGtrain = EEG;
if ~isfield(EEGtrain, 'etc') || ~isstruct(EEGtrain.etc), EEGtrain.etc = struct(); end
seg = struct('removed', false, 'n_intervals', 0, 'pct_time', 0, 'intervals', [], 'policy', policy);

if ~isfield(pre, 'initrej') || ~isfield(pre.initrej, 'badseg') || ~pre.initrej.badseg.enabled
    logmsg(logf, '[ICA-TRAIN] badseg disabled; using full data.');
    EEGtrain.etc.pp_train_seg = seg;
    return;
end

thr = pre.initrej.badseg.threshold_uv;
logmsg(logf, '[ICA-TRAIN] Detecting bad segments (thr = %.1f uV) for training copy.', thr);
badSamp   = any(abs(double(EEG.data)) > thr, 1);
intervals = mask_to_intervals(badSamp)';
seg.intervals   = intervals;
seg.n_intervals = size(intervals, 1);
seg.pct_time    = 100 * (sum(badSamp) / numel(badSamp));

if isempty(intervals)
    logmsg(logf, '[ICA-TRAIN] No bad segments detected; using full data.');
    EEGtrain.etc.pp_train_seg = seg;
    return;
end
logmsg(logf, '[ICA-TRAIN] Detected %d intervals (%.2f%% of samples).', seg.n_intervals, seg.pct_time);

switch policy
    case 'remove'
        EEGtrain = pop_select(EEGtrain, 'nopoint', intervals);
        EEGtrain = eeg_checkset(EEGtrain);
        if ~isfield(EEGtrain, 'etc') || ~isstruct(EEGtrain.etc), EEGtrain.etc = struct(); end
        seg.removed = true;
        logmsg(logf, '[ICA-TRAIN] Removed bad segments from training copy (policy = remove).');
    case 'keep'
        logmsg(logf, '[ICA-TRAIN] Keeping all segments (policy = keep).');
    otherwise
        error('st_clean_copy:Policy', 'Unknown badseg_policy "%s".', policy);
end
EEGtrain.etc.pp_train_seg = seg;
end