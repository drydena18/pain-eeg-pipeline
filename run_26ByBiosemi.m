function run_26ByBiosemi(only_sub, doICA, doPost, interactiveClean)
% run_26ByBiosemi(only_sub, doICA, doPost)
% Stage 1 (always): import BDF -> relabel -> chan locs ->
% resample/filter/notch -> epoch [-1 2] around event_marker -> baseline ->
% add CSV fields -> save to derivatives/filter/
% Stage 2 (if doICA): load from derivatives/no_ica -> runica -> save to
% derivatives/ica/ (assumes manually rejected bad chans into no_ica/
% Stage 3 (if doPost): load from derivatives/after_ica/ -> lowpass 30 ->
% save to refilter_30/ -> reref avg -> save to rerefer

if nargin < 2, doICA = false; end
if nargin < 3, doPost = false; end
if nargin < 4, interactiveClean = doICA; end

clearvars -except only_sub doICA doPost interactiveClean; clc; eeglab

P = config_paths();

% Participant table
behav_data = tdfread(P.PARTICIPANTS_TSV);
sub = behav_data.participant_id;

% Resolve only_sub to index
sub_names = strtrim(cellstr(sub));
idx = find(strcmpi(sub_names, only_sub), 1);
if isempty(idx)
    error('Participant %s not found in %s', only_sub, P.PARTICIPANTS_TSV);
end

event_marker = {'condition 54'};

% I/O Paths (stage 1)
input_path = [P.BIDS filesep];
save_path1 = [P.SESSION_MERGED filesep];
save_path2 = [P.FILTER filesep];

% Singletrial CSV read
infoT = readtable(P.CSV_SINGLETRIAL, FileType = "text");
vnames = lower(infoT.Properties.VariableNames);
idCol = infoT.Properties.VariableNames{ find(vnames=="id", 1, "first") };
trialCol = infoT.Properties.VariableNames{ find(contains(vnames,"trial") & contains(vnames,"num"), 1, "first") };
powerCol = infoT.Properties.VariableNames{ find(contains(vnames,"laser") & contains(vnames,"power"), 1, "first") };
ratingCol = infoT.Properties.VariableNames{ find(contains(vnames,"pain") & contains(vnames,"rat"), 1, "first") };

%% ==============================================
% Stage 1: Import -> Filter -> Epoch -> CSV Merge
% ===============================================

subi = idx;
input_path_1 = [input_path sub(subi, :) filesep 'eeg' filesep];
filename1 = [input_path_1, sub(subi, :) '_task-26ByBiosemi_eeg.bdf'];

EEG = pop_biosig(filename1);

% ---- modify electrode names ----
EEG = pop_select( EEG, 'channel', ...
    {'A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14','A15','A16', ...
    'A17','A18','A19','A20','A21','A22','A23','A24','A25','A26','A27','A28','A29','A30','A31','A32', ...
    'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16', ...
    'B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','B29','B30','B31','B32'});

% ---- channel location ----
EEG = pop_chanedit(EEG, 'changefield', {1 'labels' 'FP1'},'changefield',{2 'labels' 'AF7'},'changefield',{3 'labels' 'AF3'}, ...
    'changefield',{4 'labels' 'F1'},'changefield',{5 'labels' 'F3'},'changefield',{6 'labels' 'F5'}, ...
    'changefield',{7 'labels' 'F7'},'changefield',{8 'labels' 'FT7'},'changefield',{9 'labels' 'FC5'}, ...
    'changefield',{10 'labels' 'FC3'},'changefield',{11 'labels' 'FC1'},'changefield',{12 'labels' 'C1'}, ...
    'changefield',{13 'labels' 'C3'},'changefield',{14 'labels' 'C5'},'changefield',{15 'labels' 'T7'}, ...
    'changefield',{16 'labels' 'TP7'},'changefield',{17 'labels' 'CP5'},'changefield',{18 'labels' 'CP3'}, ...
    'changefield',{19 'labels' 'CP3'},'changefield',{20 'labels' 'CP1'},'changefield',{21 'labels' 'P1'}, ...
    'changefield',{22 'labels' 'P3'},'changefield',{23 'labels' 'P5'},'changefield',{24 'labels' 'P7'}, ...
    'changefield',{19 'labels' 'CP1'},'changefield',{20 'labels' 'P1'},'changefield',{21 'labels' 'P3'}, ...
    'changefield',{22 'labels' 'P5'},'changefield',{23 'labels' 'P7'},'changefield',{24 'labels' 'P9'}, ...
    'changefield',{25 'labels' 'PO7'},'changefield',{26 'labels' 'PO3'},'changefield',{27 'labels' 'O1'}, ...
    'changefield',{28 'labels' 'LZ'},'changefield',{29 'labels' 'OZ'},'changefield',{30 'labels' 'POZ'}, ...
    'changefield',{31 'labels' 'PZ'},'changefield',{32 'labels' 'CPZ'},'changefield',{33 'labels' 'FPZ'}, ...
    'changefield',{34 'labels' 'FP2'},'changefield',{35 'labels' 'AF8'},'changefield',{36 'labels' 'AF4'}, ...
    'changefield',{37 'labels' 'AFZ'},'changefield',{38 'labels' 'FZ'},'changefield',{39 'labels' 'F2'}, ...
    'changefield',{40 'labels' 'F4'},'changefield',{41 'labels' 'F6'},'changefield',{42 'labels' 'F8'}, ...
    'changefield',{43 'labels' 'FT8'},'changefield',{44 'labels' 'FC6'},'changefield',{45 'labels' 'FC4'}, ...
    'changefield',{45 'labels' 'FC4'},'changefield',{46 'labels' 'FC2'},'changefield',{47 'labels' 'FCZ'}, ...
    'changefield',{48 'labels' 'CZ'},'changefield',{49 'labels' 'C2'},'changefield',{50 'labels' 'C4'}, ...
    'changefield',{51 'labels' 'C6'},'changefield',{52 'labels' 'T8'},'changefield',{53 'labels' 'TP8'}, ...
    'changefield',{54 'labels' 'CP6'},'changefield',{55 'labels' 'CP4'},'changefield',{56 'labels' 'CP2'}, ...
    'changefield',{57 'labels' 'P2'},'changefield',{58 'labels' 'P4'},'changefield',{59 'labels' 'P6'}, ...
    'changefield',{60 'labels' 'P8'},'changefield',{61 'labels' 'P10'},'changefield',{62 'labels' 'PO8'}, ...
    'changefield',{63 'labels' 'PO4'},'changefield',{64 'labels' 'O2'}, 'lookup', P.ELP_FILE);

EEG = save_eeg_set(EEG, save_path1, [sub(subi,:) '_26ByBiosemi.set'], [strtrim(sub(subi,:)) '_merged']);

% ---- resample & filters ----
EEG = pop_resample( EEG, 1000);
EEG = pop_eegfiltnew(EEG, 'locutoff', 1); 
EEG = pop_eegfiltnew(EEG, 'hicutoff', 100); 
EEG = pop_eegfiltnew(EEG, 'locutoff', 48, 'hicutoff', 52, 'revfilt', 1); % notch filter

% ---- epoch and baseline ----
EEG = pop_epoch( EEG, event_marker, [-1, 2], 'newname', 'Merged datsets epochs', 'epochinfo', 'yes'); %epoching
EEG = pop_rmbase( EEG, [-1000 0]);

% ---- add rating and laser power from csv ----
subjid = strtrim(sub(subi,:));
csvIDs = string(infoT.(idCol));
candIDs = unique([ string(subjid), erase(string(subjid),"sub-"), extractAfter(string(subjid), "sub-") ]);

if isnumeric(infoT.(idCol))
    numID = str2double(extractAfter(string(subjid), "sub-"));
    mask = infoT.(idCol) == numID;
else
    mask = ismember(lower(strtrim(csvIDs)), lower(strtrim(candIDs)));
end

Ts = infoT(mask, :);

% sort by trial number to align with epoch order
if ~issorted(Ts.(trialCol))
    Ts = sortrows(Ts, trialCol, "ascend");
end

% find eeg.event index at latency ~0
nEp = numel(EEG.epoch);
idx0 = zeros(1,nEp);
for k = 1:nEp
    elats = EEG.epoch(k).eventlatency; if iscell(elats), elats = cell2mat(elats); end
    einds = EEG.epoch(k).event; if iscell(einds), einds = cell2mat(einds); end
    zpos = find(abs(elats) <1e-6,1, 'first'); if isempty(zpos), zpos = 1; end
    idx0(k) = einds(zpos);
end

if height(Ts) ~= nEp
    warning('CSV trials (%d) =/= epochs (%d) for %s. Assigning min.', height(Ts), nEp, subjid);
end
nWrite = min(height(Ts), nEp);

for k = 1:nWrite
    EEG.event(idx0(k)).laser_power = Ts.(powerCol)(k);
    EEG.event(idx0(k)).rating = Ts.(ratingCol)(k);
    EEG.event(idx0(k)).trial_num = Ts.(trialCol)(k);
end

    EEG = save_eeg_set(EEG, save_path2, [sub(subi,:) '_26ByBiosemi.set'], [strtrim(sub(subi,:)) '_filter']);

%% ======================
% Stage 2: ICA (Optional)
% =======================

if doICA
    file_name = [sub(subi,:) '_26ByBiosemi.set'];

    EEG = pop_loadset('filename', file_name, 'filepath', save_path2);

    if interactiveClean
        fprintf('\n*** Interactive cleaning for %s ***\n', subjid);
        fprintf('1) A viewer will open. Mark bad epochs, then close it.\n');
        fprintf('2) You will be prompted to list bad channels to REMOVE (labels, comma-separated).\n');
        fprintf(' These will be interpolated back AFTER ICA in Stage 3.\n');

        % ---- Open scroll & wait until closed ----
        pop_eegplot(EEG, 1,1,1);
        fig = findall(0, 'type', 'figure', 'tag', 'EEGPLOT');
        if isempty(fig), pause(0.1); fig = findall(0, 'type', 'figure', 'tag', 'EEGPLOT'); end
        if ~isempty(fig), waitfor(fig(1)); else, warning('EEGPLOT window not detected; continuing.'); end

        % ---- Apply manual epoch rejections ----
        EEG = eeg_rejsuperpose(EEG, 1,1,1,1,1,1,1,1);
        if isfield(EEG.reject,'rejmanual') && any(EEG.reject.rejmanual)
            EEG = pop_rejepoch(EEG, EEG.reject.rejmanual, 0);
        end

        % Ask user which channels to REMOVE before ICA (no interp yet) ----
        prompt = {'Bad channel labels to remove (comma-separated, e.g., "FP1, F8:'};
        dlgtitle = 'Remove bad channels before ICA';
        answer = inputdlg(prompt, dlgtitle, [1 80], {' '});
        badLabels = {};
        allLabs = {EEG.chanlocs.labels};
        if ~isempty(answer)
            raw = strtrim(answer{1});
            if ~isempty(raw)
                parts = regexp(raw, '\s*, \s*', 'split');
                badLabels = strtrim(parts);
                % Map labels -> indices present in current data
                allLabs = {EEG.chanlocs.labels};
                rmIdx = find(ismember(upper(allLabs), upper(badLabels)));
                if ~isempty(rmIdx)
                    % Save original chanlocs (needed for later)
                    orig_chanlocs = EEG.chanlocs;
                    % Store metadata for Stage 3
                    badinfo.badLabels = allLabs(rmIdx);
                    badinfo.subjid = subjid;
                    badinfo.timestamp = datetime("now");
                    badinfo.file = fullfile(P.NO_ICA, ['badchans_' subjid '.mat']);
                    if ~exist(P.NO_ICA, 'dir'), mkdir(P.NO_ICA); end
                    save(badinfo_file, 'badinfo', 'orig_chanlocs');

                    % Remove channels before ICA
                    EEG = pop_select(EEG, 'nochannel', rmIdx);
                else
                    warning('No matching bad channel labels found in current data; nothing removed.');
                end
            end
        end
    end

    % ---- Make sure data are double precision and consistent ----
    EEG = eeg_checkset(EEG);
    if ~isa(EEG.data, 'double'), EEG.data = double(EEG.data); end

    % Save cleaned to no_ica/
    EEG = save_eeg_set(EEG, [P.NO_ICA filesep], file_name, [strtrim(sub(subi,:)) '_clean_nointerp']);

    % ---- Run ICA; use PCA to avoid rank-mismatch warnings ----
    % A safe PCA rank is the matrix rank of the current data (channels *
    % samples)
    dat = reshape(EEG.data, EEG.nbchan, []); % channels * (time*trials)
    try
        rnk = rank(dat);
    catch
        % fallback if rank() is unavailable/slow
        rnk = min(EEG.nbchan, size(dat,2));
    end
    % Cap PCA to number of channels to be safe
    rnk = min(rnk, EEG.nbchan);

    EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', rnk, 'interrupt', 'on');

    % ---- Save ICA dataset to ica/ ----
    EEG = save_eeg_set(EEG, [P.ICA filesep], [sub(subi,:) '_26ByBiosemi.set'], [strtrim(sub(subi,:)) '_ica']);
end

%% ================================================
% Stage 3: Refilter (30Hz) + Rereference (optional)
% =================================================

if doPost
    file_name = [sub(subi,:) '_26ByBiosemi.set'];

    % Prefer dataset AFTER manual IC rejection; else fall back to ICA
    % output
    src_post = [P.AFTER_ICA filesep];
    if ~exist(fullfile(src_post, file_name), 'file')
        warning('No file in after_ica/. Falling back to ica/ for post-ICA on %s.', subjid);
        src_post = [P.ICA filesep];
    end

    EEG = pop_loadset('filename', file_name, 'filepath', src_post);

    % ---- If removed channels earlier, interpolate them back now ----
    badinfo_file = fullfile(P.NO_ICA, ['badchans_' subjid '.mat']);
    if exist(badinfo_file, 'file')
        S = load(badinfo_file); % contains badinfo, orig_chanlocs
        if isfield(S, 'orig_chanlocs') && ~isempty(S.orig_chanlocs)
            % Add back any missing channels found in orig_chanlocs
            EEG = pop_interp(EEG, S.orig_chanlocs, 'spherical');
        else
            % Fallback: If only labels present, try to match against known
            % cap
            warning('orig_chanlocs not found; unable to auto-restore montage for interpolation');
        end
    end

    % ---- Final refilter and reref ----
    EEG = pop_eegfiltnew(EEG, 'hicutoff', 30);
    EEG = save_eeg_set(EEG, [P.REFILTER_30 filesep], file_name, [strtrim(sub(subi,:)) '_refilter30']);

    EEG = pop_reref(EEG, []);
    EEG = save_eeg_set(EEG, [P.REREFER filesep], file_name, [strtrim(sub(subi,:)) '_reref']);
end

end

%% ==================
% Failsafe EEG saving
% ===================

function EEG = save_eeg_set(EEG, save_dir, file_name, setname)
    if nargin < 4 || isempty(setname)
        setname = regexprep(file_name, '\.set$', '', 'ignorecase');
    end
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG.setname = setname;
    EEG = pop_saveset(EEG, ...
        'filename', file_name, ...
        'filepath', save_dir, ...
        'check', 'on', ...
        'savemode', 'onefile');
end
