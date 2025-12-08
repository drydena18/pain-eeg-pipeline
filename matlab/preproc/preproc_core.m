function preproc_core(P, cfg, subjid, infoT, idCol, trialCol, powerCol, ratingCol, ...
    doICA, doPost, interactiveClean)
% PREPROC_CORE Single-subject preprocessing for an EEG experiment
%
% P : struct from config_paths(exp_id)
% cfg : struct from load_cfg(expXX.json)
% subjid : e.g., 'sub-01'
%
% infoT, idCol, trialCol, powerCol, ratingCol : single-trial CSV info
% doICA, doPost, interactiveClean : logical flags

    if nargin < 9 || isempty(doICA), doICA = false; end
    if nargin < 10 || isempty(doPost), doPost = false; end
    if nargin < 11 || isempty(interactiveClean), interactiveClean = doICA; end

    %% =========================================================
    % Stage 1: Import -> Relabel -> Filter -> Epoch -> CSV Merge
    % ==========================================================

    % Build BIDS input path & filename
    input_path_1 = fullfile(P.BIDS, subjid, 'eeg');
    bdf_file = fullfile(input_path_1, ...
        sprintf('%s_task-%s_eeg.bdf', subjid, cfg.experiment_id));

    fprintf('\n[Stage 1] Loading BDF for %s: \n %s\n', subjid, bdf_file);
    EEG = pop_biosig(bdf_file);

    % Select 64 Biosemi channels (A1-32 B1-32)
    EEG = pop_select( EEG, 'channel', ...
        {'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14', 'A15', 'A16', ...
        'A17', 'A18', 'A19', 'A20', 'A21', 'A22', 'A23', 'A24', 'A25', 'A26', 'A27', 'A28', 'A29', 'A30', 'A31', 'A32', ...
        'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12', 'B13', 'B14', 'B15', 'B16', ...
        'B17', 'B18', 'B19', 'B20', 'B21', 'B22', 'B23', 'B24', 'B25', 'B26', 'B27', 'B28', 'B29', 'B30', 'B31', 'B32'});

    % Channel locations / relabel 
    EEG = pop_chanedit(EEG, ...
        'changefield', {1, 'labels' 'FP1'}, 'changefield', {2, 'labels' 'AF7'}, ...
        'changefield', {3, 'labels' 'AF3'}, 'changefield', {4, 'labels' 'F1'}, ...
        'changefield', {5, 'labels' 'F3'}, 'changefield', {6, 'labels' 'F5'}, ...
        'changefield', {7, 'labels' 'F7'}, 'changefield', {8, 'labels' 'FT7'}, ...
        'changefield', {9, 'labels' 'FC5'}, 'changefield', {10, 'labels' 'FC3'}, ...
        'changefield', {11, 'labels' 'FC1'}, 'changefield', {12, 'labels' 'C1'}, ...
        'changefield', {13, 'labels' 'C3'}, 'changefield', {14, 'labels' 'C5'}, ...
        'changefield', {15, 'labels' 'T7'}, 'changefield', {16, 'labels' 'C5'}, ...
        'changefield', {17, 'labels' 'CP5'}, 'changefield', {18, 'labels' 'CP3'}, ...
        'changefield', {19, 'labels' 'CP1'}, 'changefield', {20, 'labels' 'P1'}, ...
        'changefield', {21, 'labels' 'P3'}, 'changefield', {22, 'labels' 'P5'}, ...
        'changefield', {23, 'labels' 'P7'}, 'changefield', {24, 'labels' 'P9'}, ...
        'changefield', {25, 'labels' 'PO7'}, 'changefield', {26, 'labels' 'PO3'}, ...
        'changefield', {27, 'labels' 'O1'}, 'changefield', {28, 'labels' 'Lz'}, ...
        'changefield', {29, 'labels' 'Oz'}, 'changefield', {30, 'labels' 'POz'}, ...
        'changefield', {31, 'labels' 'Pz'}, 'changefield', {32, 'labels' 'CPz'}, ...
        'changefield', {33, 'labels' 'FPz'}, 'changefield', {34, 'labels' 'FP2'}, ...
        'changefield', {35, 'labels' 'AF8'}, 'changefield', {36, 'labels' 'AF4'}, ...
        'changefield', {37, 'labels' 'AFz'}, 'changefield', {38, 'labels' 'Fz'}, ...
        'changefield', {39, 'labels' 'F2'}, 'changefield', {40, 'labels' 'F4'}, ...
        'changefield', {41, 'labels' 'F6'}, 'changefield', {42, 'labels' 'F8'}, ...
        'changefield', {43, 'labels' 'FT8'}, 'changefield', {44, 'labels' 'FC6'}, ...
        'changefield', {45, 'labels' 'FC4'}, 'changefield', {46, 'labels' 'FC2'}, ...
        'changefield', {47, 'labels' 'FCz'}, 'changefield', {48, 'labels' 'Cz'}, ...
        'changefield', {49, 'labels' 'C2'}, 'changefield', {50, 'labels' 'C4'}, ...
        'changefield', {51, 'labels' 'C6'}, 'changefield', {52, 'labels' 'T8'}, ...
        'changefield', {53, 'labels' 'TP8'}, 'changefield', {54, 'labels' 'CP6'}, ...
        'changefield', {55, 'labels' 'CP4'}, 'changefield', {56, 'labels' 'CP2'}, ...
        'changefield', {57, 'labels' 'P2'}, 'changefield', {58, 'labels' 'P4'}, ...
        'changefield', {59, 'labels' 'P6'}, 'changefield', {60, 'labels' 'P8'}, ...
        'changefield', {61, 'labels' 'P10'}, 'changefield', {62, 'labels' 'PO8'}, ...
        'changefield', {63, 'labels' 'PO4'}, 'changefield', {64, 'labels' 'O2'}, ...
        'lookup', P.ELP_FILE);

    % Base filename derived from subject + capsize
    cap_n = EEG.nbchan;
    cap_str = sprintf('%d', cap_n);
    base_name = sprintf('%s_%s', subjid, cap_str);

    % I/O Paths
    sess_path = [P.SESSION_MERGED filesep];
    filter_path = [P.FILTER filesep];
    outpath = [P.OUTPATH filesep];

    % Save merged (post-chanlocs, pre-filters)
    EEG = save_eeg_set(EEG, sess_path, [base_name '_merged.set'], [subjid '_merged']);

    % Re-reference via interpolation
    EEG = rereference_interp(EEG, 1);
    com = sprintf('Reference data');
    EEG = eeg_hist(EEG, com);
    pop_saveh(EEG.history, 'EEG_history.txt', outpath);

    % Resample & band-pass filters (config-driven)
    fprinmtf('[Stage 1] Resampling to %d Hz...\n', cfg.filters.resample_Hz);
    EEG = pop_resample(EEG, cfg.filters.resample_Hz);

    bandpass = [cfg.filters.highpass cfg.filters.lowpass];

    fprintf('[Stage 1] Low-pass filtering below %g Hz...\n', bandpass(2));
    EEG = pop_eegfilt(EEG, bandpass(2), 0, [], 0, 0, 0, 'fir1', 0);

    fprintf('[Stage 1] High-pass filtering above %g Hz...\n', bandpass(1));
    EEG = pop_eegfilt(EEG, 0, bandpass(1), 0, [], 0, 0, 0, 'fir1', 0);

    % Interactive bad-channel interpolation
    iter = 1;
    EEG.filename = sprintf('%s_%s.set', EEG.filename(1:end-4), 'initrej');
    pop_eegplot(EEG, 1, 1, 1);

    varchan = squeeze(var(EEG.data, [], 2));
    figure(30); bar(varchan);
    meanvarchan = mean(varchan);
    stdvarchan = std(varchan);
    threshdelchan = meanvarchan - (stdvarchan * 2);
    delchan = find(varchan >= threshdelchan);
    EEG.chanInterp = delchan;

    selchan = input('Channels to interpolate (labels, space- or comma-separated): ', 's');
    selchan = regexp(selchan, '(\S+)', 'tokens');
    selchan = cellfun(@(x) x{1}, selchan, 'UniformOutput', false);
    chnloc = find(ismember({EEG.chanlocs.labels}, selchan);
    if ~isempty(chnloc)
        EEG = eeg_interp(EEG, chnloc);
    end
    interpchan{iter} = chnloc; %#ok<NASGU>

    ManIsn = input(['Are you finished manually inspecting your data? ' ...
        'Press 1 to keep interpolating; Press 2 if you are finished. '], 's');
    ManIsn = str2double(ManIsn);

    while ManIsm == 1
        iter = iter + 1;
        fprintf('Continue analysis...\n');
        varchan = squeeze(var(EEG.data, [], 2));
        figure(30); bar(varchan);
        meanvarchan = mean(varchan);
        stdvarchan = std(varchan);
        threshdelchan = meanvarchan - (stdvarchan * 2);
        delchan = find(varchan >= threshdelchan);
        EEG.chanInterp = delchan;

        selchan = input('Channels to interpolate (labels, space- or comma-separated): ', 's');
        if isempty(selchan)
            fprintf('No channels selected for interpolation.\n');
        else
            selchan = regexp(selchan, '(\S+)', 'tokens');
            selchan = cellfun(@(x) x{1}, selchan, 'UniformOutput', false);
            chnloc = find(ismember({EEG.chanlocs.labels}, selchan));
            EEG = eeg_interp(EEG, chnloc);
            interpchan{iter} = chnloc; %#ok<NASGU>
        end

        ManIsn = input(['Are you finished manually inspecting your data? ' ...
            'Press 1 to keep interpolating; Press 2 if you are finished. '], 's');
        ManIsn = str2double(ManIsn);
    end

    % Epoch and Baseline (config-driven)
    event_marker = cfg.events.stim_codes;
    epoch_win = [cfg.epoching.tmin cfg.epoching.tmax];
    base_sec = cfg.epoching.baseline;
    base_ms = base_sec * 1000;

    EEG = pop_epoch(EEG, event_marker, epoch_win, ...
        'newname', 'Merged datasets epochs', ...
        'epochinfo', 'yes');
    
    EEG = pop_rmbase(EEG, base_ms);

    % Add rating and laser power from CSV
    csvIDs = string(infoT.(idCol));
    candIDs = unique([ string(subjid), ...
        erase(string(subjid), 'sub-'), ...
        extractAfter(string(subjid), 'sub-') ]);

    if isnumeric(infoT.(idCol))
        numID = str2double(extractAfter(string(subjid), 'sub-'));
        mask = infoT.(idCol) == numID;
    else
        mask = ismember(lower(strtrim(csvIDs)), lower(strtrim(candIDs)));
    end

    Ts = infoT(mask, :);

    if ~issorted(Ts.(trialCol))
        Ts = sortrows(Ts, trialCol, 'ascend');
    end

    nEp = numel(EEG.epoch);
    idx0 = zeros(1, nEp);
    for k = 1:nEp
        elats = EEG.epoch(k).eventlatency; if iscell(elats), elats = cell2mat(elats); end
        einds = EEG.epoch(k).event; if iscell(einds), einds = cell2mat(einds); end
        zpos = find(abs(elats) < 1e-6, 1, 'first'); if isempty(zpos), zpos = 1; end
        idx0(k) = einds(zpos);
    end

    if height(Ts) ~= nEp
        warning('CSV trials (%d) =/= epochs (%d) for %s. Assigning min.', height(Ts), nEp, subjid);
    end
    nWrite = min(height(Ts).nEp);

    for k = 1:nWrite
        EEG.event(idx0(k)).laser_power = Ts.(powerCol)(k);
        EEG.event(idx0(k)).rating = Ts.(ratingCol)(k);
        EEG.event(idx0(k)).trial_num = Ts.(trialCol)(k);
    end

    % Save preprocessed epochs to FILTER dir
    preproc_file = [base_name '_preproc.set'];
    EEG = save_eeg_set(EEG, filter_path, preproc_file, [subjid '_preproc']);

    %% ======================
    % Stage 2: ICA (Optional)
    % =======================

    if doICA
        fprintf('\n[Stage 2] ICA for %s\n', subjid);

        file_name = preproc_file;
        EEG = pop_loadset('filename', file_name, 'filepath', filter_path);

        log = struct();
        log.sub = subjid;
        log.when = datetime('now');
        log.stage = 'pre-ICA';

        if interactiveClean
            fprintf('\n*** Interactive cleaning for %s ***\n', subjid);
            fprintf('1) A viewer will open. Mark bad epochs, then close it.\n');
            fprintf('2) You will be prompted to list bad channels to remove (labels, comma separated).\n');
            fprintf('   These will be interpolated back AFTER ICA in Stage 3.\n');

            % View data and wait for closure
            pop_eegplot(EEG, 1, 1, 1);
            fig = findall(0, 'type', 'figure', 'tag', 'EEGPLOT');
            if isempty(fig), pause(0.1); fig = findall(0, 'type', 'figure', 'tag', 'EEGPLOT'); end
            if ~isempty(fig), waitfor(fig(1)); else, warning('EEGPLOT window not detected, continuing.'); end

            % Apply manual epoch rejections
            EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
            rejIdx = [];
            if isfield(EEG.reject, 'rejmanual') && any(EEG.reject.rejmanual)
                rejIdx = find(EEG.reject.rejmanual(:)');
                EEG = pop_rejepoch(EEG, EEG.reject.rejmanual, 0);
            end

            log.rejected_epochs = rejIdx;
            log.n_epochs_after = EEG.trials;

            % Channels to remove before ICA
            prompt = {'Bad channel labels to remove (comma-separated): '};
            dlgtitle = 'Remove bad channels before ICA';
            answer = inputdlg(prompt, dlgtitle, [1 80], {' '});
            log.bad_channels_removed = {};

            if ~isempty(answer)
                raw = strtrim(answer{1});
                if ~isempty(raw)
                    parts = regexp(raw, '\s*, \s*', 'split');
                    badLabels = strtrim(parts);
                    allLabs = {EEG.chanlocs.labels};
                    rmIdx = find(ismember(upper(allLabs), upper(badLabels)));

                    if ~isempty(rmIdx)
                        orig_chanlocs = EEG.chanlocs;
                        log.bad_channels_removed = allLabs(rmIdx);
                        if ~exist(P.LOGS, 'dir'), mkdir(P.LOGS); end
                        save(fullfile(P.LOGS, ['badchans_' subjid '.mat']), ...
                            'orig_chanlocs', 'badLabels');
                        EEG = pop_select(EEG, 'nochannel', rmIdx);
                    else
                        warning('No matching bad channel labels found; nothing removed.');
                    end
                end
            end
        end

        % Sanity
        EEG = eeg_checkset(EEG);
        if ~isa(EEG.data, 'double'), EEG.data = double(EEG.data); end

        % Save cleaned (no interp yet) to no_ica/
        EEG = save_eeg_set(EEG, [P.NO_ICA filesep], file_name, [subjid '_clean_nointerp']);

        % Run ICA with PCA dimension
        dat = reshape(EEG.data, EEG.nbchan, []); % channels x (time * trials)
        try
            rnk = rank(dat)
        catch
            rnk = min(size(dat, 1), size(dat, 2));
        end
        rnk = min(rnk, EEG.nbchan);
        log.ica_rank_used = rnk;
        lob.nbchan_before_ica = EEG.nbchan;
        log.nsamples_total = size(dat, 2);

        EEG = pop_runica(EEG, 'icatype', 'runica', 'extende', 1, ...
            'pca', rnk, 'interrupt', 'on');

        % ICLabel & Quick Reports
        try
            EEG = iclabel(EEG);
            cls = EEG.etc.ic_classification.ICLabel;
            log.iclabel_classes = cls.classes;
            log.iclabel_probs = cls.classification;

            outdir = fullfile(P.DERIV, 'ica_reports');
            if ~exist(outdir, 'dir'), mkdir(outdir); end

            EyeIdx = strcmpi(cls.classes, 'Eye');
            MusIdx = strcmpi(cls.classes, ' Muscle');
            p = cls.classification;
            cand_eye = find(p(:, EyeIdx) > 0.90);
            cand_mus = find(p(:, MusIdx) > 0.90);

            if ~isempty(cand_eye)
                f1 = figure('Visible', 'off');
                pop_topoplot(EEG, 0, cand_eye, 'Eye ICs (p > .9)', 0, 'electrodes', 'off');
                exportgraphics(f1, fullfile(outdir, [subjid '_eyeICs.png']), 'Resolution', '200');
                close(f1);
            end
            if ~isempty(cand_mus)
                f2 = figure('Visible', 'off');
                pop_topoplot(EEG, 0, cand_mus, 'Muscle ICs (p > .9)', 0, 'electrodes', 'off');
                exportgraphics(f2, fullfile(outdir, [subjid '_muscleICs.png']), 'Resolution', '200');
            end
        catch
            warning('ICLabel not available; skipping auto labels/report.');
            log.iclabel_classes = {};
            log.iclabel_probs =[];
        end

        % Save ICA dataset to ICA/
        EEG = save_eeg_set(EEG, [P.ICA filesep], file_name, [subjid '_ica']);

        % Write Stage 2 log
        if ~exist(P.LOGS, 'dir'), mkdir(P.LOGS); end
        json_file = fullfile(P.LOGS, [subjid '_preproc_log_stage2.json']);
        fid = fopen(json_file, 'w');
        if fid == -1
            error('Failed to open log file for writing: %s', json_file);
        end
        fprintf(fid, '%s', jsonencode(log));
        fclose(fid);

        summary = table(string(subjid), log.n_epochs_after, log.ica_rank_used, ...
            log.nbchan_before_ica, ...
            'VariableNames', {'sub', 'epochs_after', 'ica_rank', 'nbchan'});
        summ_path = fullfile(P.LOGS, 'summary_stage2.csv');
        if exist(summ_path, 'file')
            writetable(summary, summ_path, 'WriteMode', 'append');
        else
            writetable(summary, summ_path);
        end
    end

    %% ================================================
    % Stage 3: Refilter (30Hz) + Rereference (Optional)
    % =================================================

    if doPost
        fprintf('\n[Stage 3] Post-ICA processing for %s\n', subjid);

        file_name = preproc_file;

        % Prefer dataset AFTER manual IC rejection; else fall back to ICA
        src_post = [P.AFTER_ICA filesep];
        if ~exist(fullfile(src_post, file_name), 'file')
            warning('No file in after_ica/. Falling back to ica/ for post-ICA on %s.', subjid);
            src_post = [P.ICA filesep];
        end

        EEG = pop_loadset('filename', file_name, 'filepath', src_post);

        log = struct();
        log.sub = subjid;
        log.when = datetime('now');
        log.stage = 'post-ICA';

        % If removed channels earlier, interpolate them back now
        badfile = fullfile(P.LOGS, ['badchans_' subjid '.mat']);
        if exist(badfile, 'file')
            S = load(badfile);
            if isfield(S, 'orig_chanlocs') && ~isempty(S.orig_chanlocs)
                EEG = pop_interp(EEG, S.orig_chanlocs, 'spherical')l;
                lag.bad_channels_interpolated = S.badLabels(:)';
            else
                warning('Original chanlocs missing in %s; cannot restore montage for interpolation.', badfile);
                log.bad_channels_interpolated = {};
            end
        else
            log.bad_channels_interpolated = {};
        end

        % Estimate # of ICs removed
        try
            EEG_ica = pop_loadset('filename', file_name, 'filepath', [P.ICA filesep]);
            total = size(EEG_ica.icaweights, 1);
            kept = size(EEG.icaweights, 1);
            log.ic_total = total;
            log.ic_kept = kept;
            log.ic_removed = total - kept;
        catch
            log.ic_total = NaN;
            log.ic_kept = NaN;
            log.ic_removed = NaN;
        end
        
        % Final lowpass 30Hz and reref
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 30);
        EEG = save_eeg_set(EEG, [P.REFILTER_30 filesep], file_name, [subjid '_refilter30']);

        EEG = pop_reref(EEG, []);
        EEG = save_eeg_set(EEG, [P.REREFER filesep], file_name, [subjid '_reref']);

        % Write Stage 3 Log
        if ~exist(P.LOGS, 'dir'), mkdir(P.LOGS); end
        json_file = fullfile(P.LOGS, [subjid '_preproc_log_stage3.json']);
        fid = fopen(json_file, 'w');
        fprintf(fid, '%s', jsonencode(log));
        fclose(fid);

        summary = table(string(subjid), log.ic_total, log.ic_kept, log.ic_removed, ...
            'VariableNames', {'sub', 'ic_total', 'ic_kept', 'ic_removed'});
        summ_path = fullfile(P.LOGS, 'summary_stage3.csv');
        if exist(summ_path, 'file')
            writetable(summary, summ_path, 'WriteMode', 'append');
        else
            writetable(summary, summ_path);
        end
    end
end

%% ==========================
% Helper: Failsafe EEG saving
% ===========================
function EEG = save_eeg_set(EEG, save_dir, file_name, setname)
    if nargin < 4 || isempty(setname)
        setname = regexprep(file_name,  '\.set$', '', 'ignorecase');
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG.setname = setname;
    EEG = pop_saveset(EEG, ...
        'filename', file_name, ...
        'filepath', save_dir, ...
        'check', 'on', ...
        'savemode', 'onefile');
end