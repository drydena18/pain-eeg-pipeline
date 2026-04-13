function spectral_core(P, cfg)
% SPECTRAL_CORE  Trial-wise spectral features + interaction metrics + LEP + phase
% V 2.1.0
%
% V 2.1.0 changes vs V 2.0.1:
%   - Pain ratings are loaded from P.CORE.CSV_SINGLETRIAL (experiment-level
%     participants_singletrial CSV) via spec_load_singletrial_ratings,
%     instead of relying on EEG.epoch.pain_rating (which is typically absent)
%   - Ratings are passed explicitly to spec_compute_phase_metric (inline
%     Hilbert path) and to spec_compute_rcl_from_phase (09_hilbert phase)
%   - r_cl computation is now handled by the dedicated helper
%     spec_compute_rcl_from_phase for both phase loading paths
%
% New in V2 vs V1:
%   - Windowed pre/post-stim PSDs  (spec_compute_windowed_psds)
%   - Interaction metrics: BI_pre, LR_pre, CoG_pre, psi_cog, ΔERD
%     (spec_compute_interaction_metrics)
%   - Slow-alpha Hilbert phase per channel x trial:
%       1st: loads from pre-computed 09_hilbert stage if present
%       2nd: computes inline via spec_compute_phase_metric
%   - LEP per trial + GA waveforms + N2/P2 peak CSV  (spec_compute_lep)
%   - Subject-level summary CSV with TVI_alpha  (spec_write_subject_summary_csv)
%
% V2.0.1 change:
%   Section 7 now calls spec_compute_tvi_alpha and spec_compute_ga_rcl
%   before writing, keeping spec_write_subject_summary_csv a pure writer.
%
% Output structure per subject:
%   SPECTRAL/csv/   sub-XXX_spectral_chan_by_trial.csv
%                   sub-XXX_spectral_ga_by_trial.csv
%                   sub-XXX_subject_summary.csv
%   SPECTRAL/lep/   sub-XXX_lep_trials.mat
%                   sub-XXX_lep_ga.mat
%                   sub-XXX_lep_peaks.csv
%   SPECTRAL/figures/
%   SPECTRAL/tmp/
%   SPECTRAL/logs/

subs = cfg.exp.subjects(:);
 
inStage = "08_base";
if isfield(cfg, 'spectral') && isfield(cfg.spectral, 'input_stage') && strlength(string(cfg.spectral.input_stage)) > 0
    inStage = string(cfg.spectral.input_stage);
end
 
plotMode = "summary";
if isfield(cfg.spectral, 'qc') && isfield(cfg.spectral.qc, 'plot_mode')
    plotMode = string(cfg.spectral.qc.plot_mode);
end
 
psd     = cfg.spectral.psd;
alpha   = cfg.spectral.alpha;
windows = cfg.spectral.windows;
lepCfg  = cfg.spectral.lep;
phCfg   = cfg.spectral.phase;
 
foo     = cfg.spectral.fooof;
doFooof = isfield(foo,    'enabled') && logical(foo.enabled);
doLEP   = isfield(lepCfg, 'enabled') && logical(lepCfg.enabled);
doPhase = isfield(phCfg,  'enabled') && logical(phCfg.enabled);
 
for i = 1:numel(subs)
    subjid = subs(i);
 
    subRoot = fullfile(string(P.RUN_ROOT), sprintf('sub-%03d', subjid));
    outRoot = fullfile(string(P.SPEC_ROOT), sprintf('sub-%03d', subjid));
    outCSV  = fullfile(outRoot, 'csv');
    outFig  = fullfile(outRoot, 'figures');
    outTmp  = fullfile(outRoot, 'tmp');
    outLog  = fullfile(outRoot, 'logs');
    outLEP  = fullfile(outRoot, 'lep');
 
    spec_ensure_dir(outRoot);
    spec_ensure_dir(outCSV);
    spec_ensure_dir(outFig);
    spec_ensure_dir(outTmp);
    spec_ensure_dir(outLog);
    spec_ensure_dir(outLEP);
 
    logf = spec_open_log(outLog, subjid, 'spectral');
    cobj = onCleanup(@() spec_safe_close(logf));
 
    spec_logmsg(logf, '===== SPECTRAL V2 START sub-%03d =====', subjid);
    spec_logmsg(logf, 'Input stage: %s | Plot mode: %s', inStage, plotMode);
 
    % ---------------------------------------------------------------
    % Load 08_base
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
    inSet = char(inSet);
    [inFolder, inName, inExt] = fileparts(inSet);
    EEG = pop_loadset('filename', [inName inExt], 'filepath', inFolder);
    EEG = eeg_checkset(EEG);
 
    if EEG.trials <= 1
        spec_logmsg(logf, '[WARN] EEG not epoched (trials=%d). Skipping.', EEG.trials);
        continue;
    end
 
    chanLabels = spec_get_chanlabels(EEG);
    nChan      = EEG.nbchan;
    nTr        = EEG.trials;
 
    % ---------------------------------------------------------------
    % Load per-trial pain ratings from singletrial CSV
    % P.CORE.CSV_SINGLETRIAL = .../resource/participants_singletrial_<exp>.csv
    % These are passed to the phase computation for r_cl.
    % ---------------------------------------------------------------
    ratings = [];
    if isfield(P, 'CORE') && isfield(P.CORE, 'CSV_SINGLETRIAL')
        ratings = spec_load_singletrial_ratings(P.CORE.CSV_SINGLETRIAL, subjid, nTr, logf);
    else
        spec_logmsg(logf, '[RATINGS][WARN] P.CORE.CSV_SINGLETRIAL not defined; r_cl will be NaN.');
    end
 
    % ---------------------------------------------------------------
    % 1. Full-epoch Welch PSD  ->  existing spectral features
    % ---------------------------------------------------------------
    spec_logmsg(logf, '[PSD] Full-epoch Welch PSD (nChan=%d nTrials=%d)', nChan, nTr);
    [f, Pxx] = spec_compute_psd_trials(EEG, psd);  % [nChan x nFreq x nTr]
 
    spec_logmsg(logf, '[FEAT] Full-epoch alpha features...');
    featChan = spec_compute_alpha_features_from_psd(f, Pxx, alpha);
 
    gaPxx = squeeze(mean(Pxx, 1, 'omitnan'));   % [nFreq x nTr]
    if size(gaPxx, 1) ~= numel(f)
        gaPxx = gaPxx';
    end
    featGA = spec_compute_alpha_features_from_psd(f, reshape(gaPxx, [1 numel(f) nTr]), alpha);
    featGA = spec_squeeze_ga_features(featGA);
 
    % ---------------------------------------------------------------
    % 2. Pre/post-stim windowed PSDs  ->  interaction metrics
    % ---------------------------------------------------------------
    spec_logmsg(logf, '[WIN_PSD] pre=[%.2f %.2f]s  post=[%.2f %.2f]s', ...
        windows.pre_sec(1), windows.pre_sec(2), windows.post_sec(1), windows.post_sec(2));
 
    try
        [f2, prePxx, postPxx] = spec_compute_windowed_psds(EEG, psd, windows, logf);
    catch ME
        spec_logmsg(logf, '[WARN] Windowed PSD failed: %s', ME.message);
        f2 = []; prePxx = []; postPxx = [];
    end
 
    if ~isempty(f2)
        spec_logmsg(logf, '[INTERACT] Computing BI_pre, LR_pre, CoG_pre, DELTA_ERD...');
        try
            [intChan, intGA] = spec_compute_interaction_metrics(f2, prePxx, postPxx, alpha, logf);
            fn = fieldnames(intChan);
            for k = 1:numel(fn)
                featChan.(fn{k}) = intChan.(fn{k});
                featGA.(fn{k})   = intGA.(fn{k});
            end
        catch ME
            spec_logmsg(logf, '[WARN] Interaction metrics failed: %s', ME.message);
        end
    end
 
    % ---------------------------------------------------------------
    % 3. Slow-alpha Hilbert phase at t=0
    %    3A: prefer pre-computed 09_hilbert stage
    %    3B: compute inline via spec_compute_phase_metric
    %
    %    Pain ratings (loaded above from P.CORE.CSV_SINGLETRIAL) are
    %    passed to both paths so that GA r_cl is always populated.
    % ---------------------------------------------------------------
    phaseMat = [];
    rclTable = table();
 
    if doPhase
        hilbertDir = fullfile(subRoot, '09_hilbert');
        if exist(hilbertDir, 'dir')
            dH = dir(fullfile(hilbertDir, sprintf('*%03d*_hilbert_phase.mat', subjid)));
            if ~isempty(dH)
                [~, ix] = sort([dH.datenum], 'descend');
                hMatPath = fullfile(dH(ix(1)).folder, dH(ix(1)).name);
                try
                    hData     = load(hMatPath, 'phase_slow', 'stimOnsetIdx');
                    t0idx     = hData.stimOnsetIdx;
                    phaseFull = double(hData.phase_slow);   % [nChan x nTime x nTr]
                    phaseMat  = squeeze(phaseFull(:, t0idx, :));
                    % Guard: squeeze collapses dims when nChan=1 or nTr=1
                    if isvector(phaseMat)
                        if nChan == 1
                            phaseMat = reshape(phaseMat, 1, nTr);
                        else
                            phaseMat = reshape(phaseMat, nChan, 1);
                        end
                    end
                    spec_logmsg(logf, '[PHASE] Loaded from 09_hilbert: %s', hMatPath);
 
                    % r_cl from pre-loaded phaseMat + CSV ratings
                    if ~isempty(ratings)
                        rclTable = spec_compute_rcl_from_phase( ...
                            phaseMat, ratings, chanLabels, phCfg, logf, P, subjid);
                    else
                        spec_logmsg(logf, '[PHASE] No ratings available; r_cl skipped for 09_hilbert path.');
                    end
 
                catch ME
                    spec_logmsg(logf, '[PHASE][WARN] Load failed: %s', ME.message);
                    phaseMat = [];
                end
            end
        end
 
        if isempty(phaseMat)
            spec_logmsg(logf, '[PHASE] 09_hilbert not found; computing inline...');
            try
                % ratings passed explicitly; spec_compute_phase_metric uses
                % them instead of the EEG.epoch fallback.
                [phaseMat, rclTable] = spec_compute_phase_metric(EEG, alpha, phCfg, logf, ratings);
            catch ME
                spec_logmsg(logf, '[PHASE][WARN] Inline computation failed: %s', ME.message);
                phaseMat = [];
            end
        end
 
        if ~isempty(phaseMat) && isequal(size(phaseMat), [nChan nTr])
            featChan.phase_slow_rad = phaseMat;
            % Circular mean across channels for GA
            featGA.phase_slow_rad = angle(mean(exp(1i * phaseMat), 1, 'omitnan'))';
        else
            if ~isempty(phaseMat)
                spec_logmsg(logf, '[PHASE][WARN] phaseMat size %s != [%d %d]; discarding.', ...
                    mat2str(size(phaseMat)), nChan, nTr);
            end
        end
    end
 
    % ---------------------------------------------------------------
    % 4. FOOOF on GA PSD
    % ---------------------------------------------------------------
    fooofOut = struct();
    if doFooof
        try
            spec_logmsg(logf, '[FOOOF] Running Python FOOOF bridge...');
            fooofOut = spec_run_fooof_python(f, gaPxx, foo, outTmp, subjid, logf);
            if isfield(fooofOut, 'trials') && ~isempty(fooofOut.trials)
                fooofOut = spec_fill_fooof_alpha(fooofOut, f, gaPxx, featGA, cfg.spectral.fooof.alpha_band_hz);
            end
            spec_logmsg(logf, '[FOOOF] Done. Trials fit: %d', numel(fooofOut.trials));
        catch ME
            spec_logmsg(logf, '[WARN] FOOOF failed: %s', ME.message);
            fooofOut = struct();
        end
    else
        spec_logmsg(logf, '[FOOOF] disabled.');
    end
 
    % ---------------------------------------------------------------
    % 5. LEP per trial + GA
    % ---------------------------------------------------------------
    if doLEP
        spec_logmsg(logf, '[LEP] Extracting Laser Evoked Potentials...');
        try
            spec_compute_lep(EEG, lepCfg, outLEP, subjid, logf);
        catch ME
            spec_logmsg(logf, '[WARN] LEP computation failed: %s', ME.message);
        end
    end
 
    % ---------------------------------------------------------------
    % 6. Write CSVs
    % ---------------------------------------------------------------
    spec_logmsg(logf, '[CSV] Writing channel-by-trial and GA-by-trial CSVs...');
    outChanCSV = fullfile(outCSV, sprintf('sub-%03d_spectral_chan_by_trial.csv', subjid));
    outGaCSV   = fullfile(outCSV, sprintf('sub-%03d_spectral_ga_by_trial.csv',  subjid));
 
    spec_write_chan_trial_csv(outChanCSV, subjid, chanLabels, featChan);
    spec_write_ga_trial_csv(outGaCSV, subjid, featGA, fooofOut);
 
    % ---------------------------------------------------------------
    % 7. Subject-level summary CSV  (TVI_alpha + GA r_cl)
    %
    % spec_compute_tvi_alpha  -> tviOut    struct
    % spec_compute_ga_rcl     -> gaRclOut  struct
    % ---------------------------------------------------------------
    if isfield(featGA, 'bi_pre')
        outSumCSV = fullfile(outCSV, sprintf('sub-%03d_subject_summary.csv', subjid));
        try
            tviOut   = spec_compute_tvi_alpha(featGA.bi_pre, logf);
            gaRclOut = spec_compute_ga_rcl(rclTable, logf);
            spec_write_subject_summary_csv(outSumCSV, subjid, tviOut, gaRclOut, logf);
        catch ME
            spec_logmsg(logf, '[WARN] Subject summary CSV failed: %s', ME.message);
        end
    end
 
    % ---------------------------------------------------------------
    % 8. Trial-spectral QC plots (optional)
    % ---------------------------------------------------------------
    if isfield(cfg.spectral, 'trial_spectral') && ...
       isfield(cfg.spectral.trial_spectral, 'enabled') && ...
       logical(cfg.spectral.trial_spectral.enabled)
 
        outTrialSpec = fullfile(outRoot, 'trial_spectral');
        spec_ensure_dir(outTrialSpec);
        spec_logmsg(logf, '[TRIALSPEC] Saving pre/post-stim spectral QC plots...');
        spec_plot_trial_spectral_qc(outTrialSpec, EEG, cfg, subjid, logf);
    end
 
    % ---------------------------------------------------------------
    % 9. Summary + debug figures
    % ---------------------------------------------------------------
    try
        spec_logmsg(logf, '[FIG] Writing summary figures...');
        outSummaryFig = fullfile(outFig, sprintf('sub-%03d_spectral_summary.png', subjid));
        spec_plot_summary(outSummaryFig, f, gaPxx, featGA, fooofOut, alpha, cfg, subjid);
 
        if isfield(cfg.spectral.qc, 'save_heatmaps') && logical(cfg.spectral.qc.save_heatmaps)
            % V2: one call produces three separate files:
            %   sub-NNN_heatmap_prestim.png  (full-epoch + pre-stim metrics)
            %   sub-NNN_heatmap_poststim.png (ERD metrics)
            %   sub-NNN_heatmap_phase.png    (phase heatmap + rose)
            spec_plot_heatmap_panel(outFig, featChan, chanLabels, subjid);
        end
 
        if isfield(featGA, 'bi_pre')
            outInteractFig = fullfile(outFig, sprintf('sub-%03d_interaction_summary.png', subjid));
            spec_plot_interaction_summary(outInteractFig, featChan, featGA, chanLabels, subjid, logf);
        end
 
        if plotMode == "debug" || plotMode == "exhaustive"
            spec_logmsg(logf, '[FIG] Debug mode: per-trial plots...');
            spec_plot_debug_trials(outFig, f, Pxx, featChan, chanLabels, alpha, cfg, subjid, plotMode);
        end
    catch ME
        spec_logmsg(logf, '[WARN] Plotting failed: %s', ME.message);
    end
 
    spec_logmsg(logf, '===== SPECTRAL V2 DONE sub-%03d =====', subjid);
end
end