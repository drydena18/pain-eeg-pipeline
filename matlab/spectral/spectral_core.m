function spectral_core(P, cfg)
% SPECTRAL_CORE  Channel-space alpha features, computed exactly as in source space
% V 3.1.0
%
% V3.1.0: whole window = features.windows.whole_sec; slow_alpha_frac dropped;
% delta_erd -> erd_asym (see helpers).
%
% V3.0.0 changes vs V2.2.0 (channel/source parity; Python twin is
% python/source/source_core.py V4.0.0):
%   - Band power: filter-Hilbert (spec_band_power_tcs) with the shared FIR,
%     replacing Welch-trapezoid band power on 0.70 s segments, whose ~1 Hz
%     bins left part of 8-12 Hz outside both sub-bands.
%   - CoG / FOOOF PSD: shared Welch (spec_welch_psd), periodic Hamming,
%     zero-padded to 0.25 Hz, each window segmented on its own.
%   - All feature parameters from cfg.spectral.feat (= cfg.features).
%   - Spatial GA per trial is now consistent across windows (ratio of mean
%     powers everywhere; V2.2.0 used mean-of-ratios for pre/post) and its
%     p5_flag comes from the GA unit's own powers (V2.2.0 OR-ed ~64 channels,
%     flagging nearly every trial).
%   - ERD guard eps0 / p5 thresholds per channel, data-relative (no 45-55 Hz
%     quiet band, which the 40 Hz low-pass removes).
%   - Phase: inline with the shared FIR at t = 0 and at features.phase.post_ref_sec,
%     columns slow_phase / sin_phase / cos_phase (+ _post). The 09_hilbert
%     stage output is no longer read (its pop_eegfiltnew filter differs from
%     source space). r_cl against ratings is unchanged.
%   - LEP features (shared windows, source names) are added to the
%     trial x channel CSV.
%   - NEW sub-XXX_spectral_chan_summary.csv: subject x channel, the
%     channel-space counterpart of sub-XXX_source_ga.csv (TVI per channel,
%     ITC, FOOOF on the trial-mean PSD).
%   - FOOOF fits per channel on the trial-mean whole-epoch PSD via the shared
%     Python implementation; per-trial channel-GA fits are dropped.
%
% Output structure per subject (P.SPEC_ROOT/sub-XXX/):
%   csv/  sub-XXX_spectral_chan_by_trial.csv   trial x channel  (read by R merge)
%         sub-XXX_spectral_ga_by_trial.csv     trial x channel-mean
%         sub-XXX_spectral_chan_summary.csv    subject x channel
%         sub-XXX_subject_summary.csv          TVI of channel-mean + r_cl
%   lep/  sub-XXX_lep_trials.mat, sub-XXX_lep_ga.mat, sub-XXX_lep_peaks.csv
%   figures/ tmp/ logs/

subs = cfg.exp.subjects(:);

inStage = "08_base";
if isfield(cfg.spectral, 'input_stage') && strlength(string(cfg.spectral.input_stage)) > 0
    inStage = string(cfg.spectral.input_stage);
end

plotMode = string(cfg.spectral.qc.plot_mode);
F       = cfg.spectral.feat;
lepCfg  = cfg.spectral.lep;
phCfg   = cfg.spectral.phase;
doFooof = logical(F.fooof.enabled);
doLEP   = logical(lepCfg.enabled);
doPhase = logical(phCfg.enabled);

for i = 1:numel(subs)
    subjid = subs(i);

    subRoot = fullfile(string(P.RUN_ROOT), sprintf('sub-%03d', subjid));
    outRoot = fullfile(string(P.SPEC_ROOT), sprintf('sub-%03d', subjid));
    outCSV  = fullfile(outRoot, 'csv');
    outFig  = fullfile(outRoot, 'figures');
    outTmp  = fullfile(outRoot, 'tmp');
    outLog  = fullfile(outRoot, 'logs');
    outLEP  = fullfile(outRoot, 'lep');
    cellfun(@spec_ensure_dir, {outRoot, outCSV, outFig, outTmp, outLog, outLEP});

    logf = spec_open_log(outLog, subjid, 'spectral');
    cobj = onCleanup(@() spec_safe_close(logf)); %#ok<NASGU>

    spec_logmsg(logf, '===== SPECTRAL V3 START sub-%03d =====', subjid);
    spec_logmsg(logf, 'Input stage: %s | Plot mode: %s', inStage, plotMode);

    % ---------------------------------------------------------------
    % Load epoched data
    % ---------------------------------------------------------------
    inDir = fullfile(subRoot, char(inStage));
    if ~exist(inDir, 'dir')
        spec_logmsg(logf, '[WARN] Missing input dir: %s (skipping)', inDir);
        continue;
    end
    inSet = spec_find_latest_set(inDir, cfg.exp.out_prefix, subjid);
    if strlength(inSet) == 0
        spec_logmsg(logf, '[WARN] No .set found in %s (skipping)', inDir);
        continue;
    end
    spec_logmsg(logf, '[LOAD] %s', inSet);
    [inFolder, inName, inExt] = fileparts(char(inSet));
    EEG = pop_loadset('filename', [inName inExt], 'filepath', inFolder);
    EEG = eeg_checkset(EEG);
    if EEG.trials <= 1
        spec_logmsg(logf, '[WARN] EEG not epoched (trials=%d). Skipping.', EEG.trials);
        continue;
    end

    chanLabels = spec_get_chanlabels(EEG);
    nTr        = EEG.trials;
    fs         = EEG.srate;
    timesSec   = double(EEG.times(:))' / 1000;

    ratings = [];
    if isfield(P, 'CORE') && isfield(P.CORE, 'CSV_SINGLETRIAL')
        ratings = spec_load_singletrial_ratings(P.CORE.CSV_SINGLETRIAL, subjid, nTr, logf);
    end

    % ---------------------------------------------------------------
    % 1. All features (shared estimators)
    % ---------------------------------------------------------------
    spec_logmsg(logf, '[FEAT] fs=%g Hz, %d chans x %d trials; filter-Hilbert trans_bw=%.2f Hz', ...
        fs, EEG.nbchan, nTr, F.band_power.trans_bw_hz);
    try
        [featChan, featGA, summ, aux] = spec_compute_subject_features(EEG.data, timesSec, fs, F);
    catch ME
        spec_logmsg(logf, '[ERROR] Feature computation failed: %s', getReport(ME, 'extended', 'hyperlinks', 'off'));
        continue;
    end
    spec_logmsg(logf, '[FEAT] p5_flag: %d / %d trial x channel cells flagged', ...
        sum(featChan.p5_flag(:)), numel(featChan.p5_flag));

    % ---------------------------------------------------------------
    % 2. FOOOF: subject x channel trial-mean whole-epoch PSD
    % ---------------------------------------------------------------
    if doFooof
        try
            fo = spec_run_fooof_python(aux.f, aux.wholePsdTrialMean, F, cfg.spectral.fooof, outTmp, subjid, logf);
            summ = spec_merge_structs(summ, fo);
        catch ME
            spec_logmsg(logf, '[WARN] FOOOF failed: %s', ME.message);
        end
    end

    % ---------------------------------------------------------------
    % 3. r_cl (MATLAB-only extra): phase at t = 0 vs single-trial ratings
    % ---------------------------------------------------------------
    rclTable = table();
    if doPhase && ~isempty(ratings)
        try
            rclTable = spec_compute_rcl_from_phase(aux.phaseT0, ratings, chanLabels, phCfg, logf, P, subjid);
        catch ME
            spec_logmsg(logf, '[PHASE][WARN] r_cl failed: %s', ME.message);
        end
    end

    % ---------------------------------------------------------------
    % 4. Write CSVs
    % ---------------------------------------------------------------
    spec_write_chan_trial_csv(fullfile(outCSV, sprintf('sub-%03d_spectral_chan_by_trial.csv', subjid)), ...
        subjid, chanLabels, featChan);
    spec_write_ga_trial_csv(fullfile(outCSV, sprintf('sub-%03d_spectral_ga_by_trial.csv', subjid)), ...
        subjid, featGA, struct());
    spec_write_chan_summary_csv(fullfile(outCSV, sprintf('sub-%03d_spectral_chan_summary.csv', subjid)), ...
        subjid, chanLabels, summ);

    try
        tviOut   = spec_compute_tvi_alpha(aux.gaPreSfBalance, logf);
        gaRclOut = spec_compute_ga_rcl(rclTable, logf);
        spec_write_subject_summary_csv(fullfile(outCSV, sprintf('sub-%03d_subject_summary.csv', subjid)), ...
            subjid, tviOut, gaRclOut, logf);
    catch ME
        spec_logmsg(logf, '[WARN] Subject summary CSV failed: %s', ME.message);
    end

    % ---------------------------------------------------------------
    % 5. LEP waveforms + peak table
    % ---------------------------------------------------------------
    if doLEP
        try
            spec_compute_lep(EEG, timesSec, F, lepCfg, outLEP, subjid, logf);
        catch ME
            spec_logmsg(logf, '[WARN] LEP save failed: %s', ME.message);
        end
    end

    % ---------------------------------------------------------------
    % 6. QC figures (best effort; plotting helpers predate V3 names)
    % ---------------------------------------------------------------
    if isfield(cfg.spectral, 'trial_spectral') && logical(cfg.spectral.trial_spectral.enabled)
        try
            outTrialSpec = fullfile(outRoot, 'trial_spectral');
            spec_ensure_dir(outTrialSpec);
            spec_plot_trial_spectral_qc(outTrialSpec, EEG, cfg, subjid, logf);
        catch ME
            spec_logmsg(logf, '[WARN] Trial-spectral QC failed: %s', ME.message);
        end
    end
    try
        spec_plot_summary(fullfile(outFig, sprintf('sub-%03d_spectral_summary.png', subjid)), ...
            aux.f, aux.gaWholePsd, featGA, struct(), cfg.spectral.alpha, cfg, subjid);
        if logical(cfg.spectral.qc.save_heatmaps)
            spec_plot_heatmap_panel(outFig, featChan, chanLabels, subjid);
        end
        spec_plot_interaction_summary(fullfile(outFig, sprintf('sub-%03d_interaction_summary.png', subjid)), ...
            featChan, featGA, chanLabels, subjid, logf);
        if plotMode == "debug" || plotMode == "exhaustive"
            spec_plot_debug_trials(outFig, aux.f, aux.wholePsd, featChan, chanLabels, cfg.spectral.alpha, cfg, subjid, plotMode);
        end
    catch ME
        spec_logmsg(logf, '[WARN] Plotting failed (cosmetic): %s', ME.message);
    end

    spec_logmsg(logf, '===== SPECTRAL V3 DONE sub-%03d =====', subjid);
end
end
