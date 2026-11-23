function spectral_core(P, cfg)
% SPECTRAL_CORE  Channel-space alpha features, computed exactly as in source space
% V 3.3.0
%
% V3.3.0: cfg.spectral.units selects what is computed and written:
%   "all"        electrodes + clusters (default; identical to V3.2.0)
%   "electrodes" electrodes only (no cluster rows)
%   "clusters"   clusters only. Only the cluster member electrodes (7 of 64
%                with the default S1 + SM_furman clusters) go through feature
%                extraction, and the per-electrode extras are skipped: channel
%                FOOOF, the whole-head GA CSV and subject summary, r_cl,
%                LEP waveform files and QC figures. The two CSVs the R merge
%                reads (chan_by_trial, chan_summary) then hold cluster rows
%                only. Cluster values are identical to an "all" run (tested:
%                tests/test_spectral_clusters.m).
%   Switching mode overwrites the per-subject CSVs: after a "clusters" run
%   the merged database has no electrode rows for those subjects.
%
% V3.2.0: electrode clusters (cfg.spectral.clusters, see spectral_default).
%   Each cluster (e.g. S1_contra / S1_ipsi relative to the stimulated hand)
%   is computed as one unit by spec_cluster_features - power pooled across
%   its electrodes, waveform features from the cluster-mean signal - and
%   appended as extra rows (chan_label = cluster name) to
%   sub-XXX_spectral_chan_by_trial.csv and sub-XXX_spectral_chan_summary.csv,
%   so the R merge ingests them unchanged; select them in the GAMM stage with
%   gamm_cfg$channel_units. Channel rows are byte-for-byte unchanged.
%   New: csv/sub-XXX_cluster_qc.csv (electrodes found / missing per cluster).
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
%   csv/  sub-XXX_spectral_chan_by_trial.csv   trial x unit  (read by R merge)
%         sub-XXX_spectral_chan_summary.csv    subject x unit (read by R merge)
%         sub-XXX_cluster_qc.csv               electrodes per cluster
%         sub-XXX_spectral_ga_by_trial.csv     trial x channel-mean   (not in "clusters")
%         sub-XXX_subject_summary.csv          TVI of channel-mean + r_cl (not in "clusters")
%   lep/  sub-XXX_lep_trials.mat, sub-XXX_lep_ga.mat, sub-XXX_lep_peaks.csv (not in "clusters")
%   figures/ tmp/ logs/
%   unit = electrode and/or cluster, per spectral.units.

subs = cfg.exp.subjects(:);

inStage = "08_base";
if isfield(cfg.spectral, 'input_stage') && strlength(string(cfg.spectral.input_stage)) > 0
    inStage = string(cfg.spectral.input_stage);
end

plotMode = string(cfg.spectral.qc.plot_mode);
F       = cfg.spectral.feat;
lepCfg  = cfg.spectral.lep;
phCfg   = cfg.spectral.phase;
clCfg   = cfg.spectral.clusters;
units   = string(cfg.spectral.units);       % "all" | "electrodes" | "clusters" (spectral_default)

clusOnly  = units == "clusters";
doClus    = logical(clCfg.enabled) && units ~= "electrodes";
doElec    = ~clusOnly;                     % electrode rows + per-electrode extras
doFooof   = logical(F.fooof.enabled);
doLEP     = logical(lepCfg.enabled) && doElec;
doPhase   = logical(phCfg.enabled);

for i = 1:numel(subs)
    subjid = subs(i);
    tSub = tic;

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

    spec_logmsg(logf, '===== SPECTRAL V3.3 START sub-%03d =====', subjid);
    spec_logmsg(logf, 'Input stage: %s | Units: %s | Plot mode: %s', inStage, units, plotMode);

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

    % ---------------------------------------------------------------
    % 0. Clusters: resolve against the full montage, write QC
    % ---------------------------------------------------------------
    C = [];
    if doClus
        try
            C = spec_resolve_clusters(clCfg, chanLabels);
            for c = 1:numel(C)
                spec_logmsg(logf, '[CLUS] %-12s found %d/%d: %s%s', C(c).name, numel(C(c).idx), ...
                    numel(C(c).requested), strjoin(C(c).found, ' '), ...
                    local_missing_txt(C(c).missing, C(c).ok, clCfg.min_channels));
            end
            spec_write_cluster_qc_csv(fullfile(outCSV, sprintf('sub-%03d_cluster_qc.csv', subjid)), ...
                subjid, C, clCfg.stim_side, clCfg.min_channels);
        catch ME
            spec_logmsg(logf, '[CLUS][WARN] Could not resolve clusters: %s', ME.message);
            C = [];
        end
    end
    if clusOnly && (isempty(C) || ~any([C.ok]))
        spec_logmsg(logf, '[CLUS][WARN] units = "clusters" but no cluster has >= %d electrodes. Skipping subject.', ...
            clCfg.min_channels);
        continue;
    end

    % Electrodes that go through feature extraction
    if clusOnly
        [members, Cf] = spec_cluster_subset(C, EEG.nbchan);
        Xin    = EEG.data(members, :, :);
        labIn  = chanLabels(members);
        spec_logmsg(logf, '[CLUS] units = "clusters": extracting %d of %d electrodes (%s)', ...
            numel(members), EEG.nbchan, strjoin(cellstr(labIn(:)'), ' '));
    else
        Cf     = C;
        Xin    = EEG.data;
        labIn  = chanLabels;
    end

    ratings = [];
    if doElec && isfield(P, 'CORE') && isfield(P.CORE, 'CSV_SINGLETRIAL')
        ratings = spec_load_singletrial_ratings(P.CORE.CSV_SINGLETRIAL, subjid, nTr, logf);
    end

    % ---------------------------------------------------------------
    % 1. All features (shared estimators)
    % ---------------------------------------------------------------
    spec_logmsg(logf, '[FEAT] fs=%g Hz, %d chans x %d trials; filter-Hilbert trans_bw=%.2f Hz', ...
        fs, size(Xin, 1), nTr, F.band_power.trans_bw_hz);
    try
        [featChan, featGA, summ, aux] = spec_compute_subject_features(Xin, timesSec, fs, F);
    catch ME
        spec_logmsg(logf, '[ERROR] Feature computation failed: %s', getReport(ME, 'extended', 'hyperlinks', 'off'));
        continue;
    end
    spec_logmsg(logf, '[FEAT] p5_flag: %d / %d trial x channel cells flagged', ...
        sum(featChan.p5_flag(:)), numel(featChan.p5_flag));

    % ---------------------------------------------------------------
    % 2. FOOOF: subject x channel trial-mean whole-epoch PSD (electrodes)
    % ---------------------------------------------------------------
    if doFooof && doElec
        try
            fo = spec_run_fooof_python(aux.f, aux.wholePsdTrialMean, F, cfg.spectral.fooof, outTmp, subjid, logf);
            summ = spec_merge_structs(summ, fo);
        catch ME
            spec_logmsg(logf, '[WARN] FOOOF failed: %s', ME.message);
        end
    end

    % ---------------------------------------------------------------
    % 2b. Cluster features (+ cluster FOOOF)
    % ---------------------------------------------------------------
    featC = struct(); summC = struct(); clusNames = {};
    if doClus && ~isempty(Cf)
        try
            [featC, summC, auxC] = spec_cluster_features(Xin, timesSec, fs, F, aux, Cf);
            clusNames = auxC.names;

            if doFooof && any(auxC.ok)
                try
                    foC = spec_run_fooof_python(aux.f, auxC.wholePsdTrialMean(auxC.ok, :), F, ...
                        cfg.spectral.fooof, fullfile(outTmp, 'clusters'), subjid, logf);
                    fn = fieldnames(foC);
                    for k = 1:numel(fn)
                        v = nan(numel(Cf), 1);
                        v(auxC.ok) = foC.(fn{k});
                        summC.(fn{k}) = v;
                    end
                catch ME
                    spec_logmsg(logf, '[CLUS][WARN] Cluster FOOOF failed: %s', ME.message);
                end
            end
        catch ME
            spec_logmsg(logf, '[CLUS][WARN] Cluster features failed: %s', ...
                getReport(ME, 'extended', 'hyperlinks', 'off'));
            featC = struct(); summC = struct(); clusNames = {};
        end
    end
    if clusOnly && isempty(clusNames)
        spec_logmsg(logf, '[CLUS][WARN] units = "clusters" and no cluster features computed. Nothing written.');
        continue;
    end

    % ---------------------------------------------------------------
    % 3. r_cl (MATLAB-only extra, per electrode): phase at t = 0 vs ratings
    % ---------------------------------------------------------------
    rclTable = table();
    if doElec && doPhase && ~isempty(ratings)
        try
            rclTable = spec_compute_rcl_from_phase(aux.phaseT0, ratings, chanLabels, phCfg, logf, P, subjid);
        catch ME
            spec_logmsg(logf, '[PHASE][WARN] r_cl failed: %s', ME.message);
        end
    end

    % ---------------------------------------------------------------
    % 4. Write CSVs
    % ---------------------------------------------------------------
    if clusOnly
        trialOut = featC;  summOut = summC;  labOut = clusNames(:);
    elseif ~isempty(clusNames)
        [trialOut, labOut] = spec_append_units(featChan, chanLabels, featC, clusNames);
        summOut            = spec_append_units(summ,     chanLabels, summC, clusNames);
    else
        trialOut = featChan; summOut = summ; labOut = chanLabels;
    end
    spec_write_chan_trial_csv(fullfile(outCSV, sprintf('sub-%03d_spectral_chan_by_trial.csv', subjid)), ...
        subjid, labOut, trialOut);
    spec_write_chan_summary_csv(fullfile(outCSV, sprintf('sub-%03d_spectral_chan_summary.csv', subjid)), ...
        subjid, labOut, summOut);
    spec_logmsg(logf, '[WRITE] %d unit(s) per trial: %d electrode(s), %d cluster(s)', ...
        numel(labOut), numel(labOut) - numel(clusNames), numel(clusNames));

    if doElec
        spec_write_ga_trial_csv(fullfile(outCSV, sprintf('sub-%03d_spectral_ga_by_trial.csv', subjid)), ...
            subjid, featGA, struct());
        try
            tviOut   = spec_compute_tvi_alpha(aux.gaPreSfBalance, logf);
            gaRclOut = spec_compute_ga_rcl(rclTable, logf);
            spec_write_subject_summary_csv(fullfile(outCSV, sprintf('sub-%03d_subject_summary.csv', subjid)), ...
                subjid, tviOut, gaRclOut, logf);
        catch ME
            spec_logmsg(logf, '[WARN] Subject summary CSV failed: %s', ME.message);
        end
    end

    % ---------------------------------------------------------------
    % 5. LEP waveforms + peak table (electrodes)
    % ---------------------------------------------------------------
    if doLEP
        try
            spec_compute_lep(EEG, timesSec, F, lepCfg, outLEP, subjid, logf);
        catch ME
            spec_logmsg(logf, '[WARN] LEP save failed: %s', ME.message);
        end
    end

    % ---------------------------------------------------------------
    % 6. QC figures (electrodes; best effort, plotting helpers predate V3 names)
    % ---------------------------------------------------------------
    if doElec
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
    end

    spec_logmsg(logf, '===== SPECTRAL V3.3 DONE sub-%03d (%.1f s) =====', subjid, toc(tSub));
end
end

% -------------------------------------------------------------------------
function txt = local_missing_txt(missing, ok, minCh)
txt = '';
if ~isempty(missing)
    txt = sprintf(' | MISSING: %s', strjoin(missing, ' '));
end
if ~ok
    txt = sprintf('%s | fewer than %d electrodes -> NaN rows', txt, minCh);
end
end