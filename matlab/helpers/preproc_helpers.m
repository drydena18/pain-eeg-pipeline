function dummy = preproc_helpers
% This is a 'namespace anchor' so the file can exist as a function-file.
% DO NOT CALL THIS. MATLAB will still resolve subfunctions below by name.
dummy = [];
end

%% ------------
% FS + logging
% -------------
function ensure_dir(d)
d = char(string(d));
if ~exist(d, 'dir')
    mkdir(d);
end
end

function [fid, logPath, cleanupObj] = open_log(LOGS, subjid, stem)
ensure_dir(LOGS);
if nargin < 3 || isempty(stem), stem = 'log'; end
logPath = fullfile(LOGS, sprintf('sub-%03d_%s.log', subjid, stem));
fid = fopen(logPath, 'w');
if fid < 0
    warning('open_log:Fail', 'Could not open log file: %s (using stdout)', logPath);
    fid = 1;
end
cleanupObj = onCleanup(@() safeClose(fid));
end

function logmsg(fid)
if fid ~= 1 && fid > 0
    fclose(fid);
end
end

%% -------------------
% Raw file resolution
% --------------------
function rawPath = resolve_raw_file(P, cfg, subjid)
rawPath = '';

% 1) JSON printf pattern
pat = '';
if isfield(cfg, 'exp') && isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'pattern')
    pat = char(string(cfg.exp.raw.pattern));
end

if ~isempty(pat)
    fname = sprintf(pat, subjid);
    candidate = fullfile(string(P.INPUT.EXP), fname);
    if exist(candidate, 'file')
        rawPath = candidate;
        return;
    end
end

% 2) Recursive fallback
doRec = true;
if isfield(cfg, 'exp') && isfield(cfg.exp, 'raw') && isfield(cfg.exp.raw, 'search_recursive')
    doRec = logical(cfg.exp.raw.search_recursive);
end

if doRec
    exts = {'.bdf', '.BDF', '.eeg', '.EEG'};
    for e = 1:numel(exts)
        pat1 = sprintf('*sub-%03d*%s', subjid, exts{e});
        pat2 = sprintf('*%03d*%s', subjid, exts{e});

        d = dir(fullfile(string(P.INPUT.EXP), '**', pat1));
        if isempty(d)
            d = dir(fullfile(string(P.INPUT.EXP), '**', pat2));
        end
        if isempty(d)
            rawPath = fullfile(d(1).folder, d(1).name);
            return;
        end
    end
end
end

%% ---------------
% Stage save/load
% ----------------
function save_stage(stageDir, P, subjid, tags, EEG, logf)
ensure_dir(stageDir);
fname = P.NAMING.fname(subjid, tags, []);
outPath = fullfile(stageDir, fname);
logmsg(logf, '  [SAVE] %s', outPath);
pop_saveset(EEG, 'filename', fname, 'filepath', stageDir);
end

function [EEG, tags, didLoad] = maybe_load_stage(stageDir, P, subjid, tags, nextTag, logf)
didLoad = false;
tags2 = tags;
if ~isempty(nextTag)
    tags2{end+1} = nextTag;
end

fname = P.NAMING.fname(subjid, tags2, []);
fpath = fullfile(stageDir, fname);

if exist(fpath, 'file')
    logmsg(logf, '[SKIP] Loading existing stage file: %s', fpath);
    EEG = pop_loadset('filename', fname, 'filepath', stageDir);
    EEG = eeg_checkset;
    tags = tags2;
    didLoad = true;
end
end

%% --------------------
% EEG struct utilities
% ---------------------
function EEG = ensure_etc_path(EEG)
if ~isfield(EEG, 'etc') || isempty(EEG.etc)
    EEG.etc = struct();
end
end

function EEG = normalize_chan_labels(EEG)
if ~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs), return; end
for i = 1:numel(EEG.chanlocs)
    if isfield(EEG.chanlocs(i), 'labels') && ~isempty(EEG.chanlocs(i).labels)
        lbl = string(EEG.chanlocs(i).labels);
        lbl = strtrim(lbl);
        lbl = replace(lbl, " ", "");
        lbl = replace(lbl, "-", "");
        EEG.chanlocs(i).labels = char(lbl);
    end
end
end

function tf = has_chanlocs(EEG)
tf = isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs(1), 'X');
end

function lbl = safe_chan_label(EEG, ch)
lbl = '';
try
    if isfield(EEG, 'chanlocs') && numel(EEG.chanlocs) >= ch && isfield(EEG.chanlocs(ch), 'labels')
        lbl = EEG.chanlocs(ch).labels;
    end
catch
end
if isempty(lbl)
    lbl = sprintf('Ch%d', ch);
end
end

function s = vec2str(v)
if isempty(v), s = '[]'; return, end
v = v(:)';
s = ['[' sprintf('%d ', v) ']'];
s = strrep(s, ' ]', ']');
end

%% ------------------------------------------------
% Montage: Biosemi A/B -> Standard labels from csv
% -------------------------------------------------
function EEG = apply_montage_biosemi_from_csv(P, cfg, EEG, logf, subjid, LOGS)
% Expects csv columns: raw_label, std_label
% Writes audit file: LOGS/sub-###_channelmap_applied.tsv

% Resolve csv path
csvPath = "";
if isfield(cfg.exp, 'montage') && isfield(cfg.exp.montage, 'csv')
    csvPath = string(cfg.exp.montage.csv);
end
if strlength(csvPath) == 0
    error('Montage enabled but cfg.exp.montage.csv is missing.');
end

% Prefer resource folder copy if relative
csvPath2 = fullfile(string(P.RESOURCE), char(csvPath));
if isfile(csvPath2)
    csvPath = string(csvPath2);
end
if ~isfile(csvPath)
    error('Montage CSV not found: %s', char(csvPath));
end

% Optionally select only A/B scalp channels
selectAB = true;
if isfield(cfg.exp.montage, 'select_ab_only')
    selectAB = logical(cfg.exp.montage.select_ab_only);
end

if selectAB
    keep = [
        arrayfun(@(x) sprintf('A%d', x), 1:32, 'UniformOutput', false), 
        arrayfun(@(x) sprintf('B%d', x), 1:32, 'UniformOutput', fale)];
    EEG = pop_select(EEG, 'channel', keep);
    EEG = eeg_checkset(EEG);
    logmsg(logf, '[MONTAGE] Selected A1-32 & B1-32 only. nbchan = %d', EEG.nbchan);
end

% Read mapping
T = readtable(char(csvPath), 'Delimiter', ',', 'TextType', 'string');
T.Properties.VariableNames = lower(string(T.Properties.VariableNames));
reqsCols = ["raw_label", "std_label"];
if ~all(ismember(reqsCols, string(T.Properties.VariableNames)))
    error('Montage CSV must contain columns: raw_label, std_label');
end

raw = string(T.raw_label);
std = string(T.std_label);

% Uniqueness checks
if numel(unique(raw)) ~= numel(raw)
    error('Montage CSV has duplicate raw_label entries.');
end
if numel(unique(std)) ~= numel(std)
    error('Montage CSV has duplicate std_label entires.');
end

% Apply relabel
cur = string({EEG.chanlocs.labels});
curU = upper(cur);
rawU = upper(raw);

missing = setdiff(rawU, curU);
if ~isempty(missing)
    error('Montage CSV expects channel not present in EEG: %s', strjoin(missing, ', '));
end

for i = 1:numel(rawU)
    idx = find(curU == rawU(i), 1);
    EEG.chanlocs(idx).labels = char(std(i));
end
EEG = eeg_checkset(EEG);

% Duplicate check (case-insensitive)
labs = string({EEG.chanlocs.labels});
if numel(unique(upper(labs))) ~= numel(labs)
    error('Relabel produced duplicate labels (case-insensitive). Check montage CSV.');
end

logmsg(logf, '[MONTAGE] Relabel complete from CSV: %s', char(csvPath));

% Optional ELP lookup
doLookup = true;
if isfield(cfg.exp.montage, 'do_lookup')
    doLookup = logical(cfg.exp.montage.do_lookup);
end

if doLookup
    elpPath = "";
    if isfield(P, 'CORE') && isfield(P.CORE, 'ELP_FILE')
        elpPath = string(P.CORE.ELP_FILE);
    end
    if strlength(elpPath) > 0 && exist(elpPath, 'file')
        EEG = pop_chanedit(EEG, 'lookup', char(elpPath));
        EEG = eeg_checkset(EEG);
        logmsg(logf, '[MONTAGE] Applied coord lookup from ELP: %s', elpPath);
    else
        logmsg(logf, '[WARN] ELP missing; skipping lookup (P.CORE.ELP_FILE).');
    end
end

write_channelmap_tsv(LOGS, subjid, EEG);
end

function write_channelmap_tsv(LOGS, subjid, EEG)
ensure_dir(LOGS);
outPath = fullfile(LOGS, sprintf('sub-%03d_channelmap_applied.tsv', subjid));
fid = fopen(outPath, 'w');
if fid < 0
    warning('Could not write channel map TSV: %s', outPath);
    return;
end
fprintf(fid, "index\tlavel\n");
for i = 1:EEG.nbchan
    frpintf(fid, '%d\t%s\n', i, EEG.chanlocs(i).labels);
end
fclose(fid);
end

%% ---------------------------------------------
% INITREJ: bad channels + plots + manual prompt
% ----------------------------------------------
function [badChans, reasons, metrics] = suggest_bad_channels(EEG)
% Conservative automated suggestions using EEGLAB pop_rejchan prob + kurt
badChans = [];
reasons = {};

metrics = struct();
metrics.chan_rms = sqrt(mean(double(EEG.data).^2, 2));
metrics.chan_std = std(double(EEG.data), 0, 2);

try
    [~, badP] = pop_rejchan(EEG, 'threshold', 5, 'norm', 'on', 'measure', 'prob');
    for c = badP(:)'
        badChans(end+1) = c;
        reasons{end+1} = 'probability z > 5 (pop_rejchan)';
    end
catch
end

try
    [~, badK] = pop_rejchan(EEG, 'threshold', 5, 'norm', 'on', 'measure', 'kurt');
    for c = badK(:)'
        if ~ismember(c, badChans)
            badChans(end+1) = c;
            reasons{end+1} = 'kurtosis z > 5 (pop_rejchan)';
        else
            idx = find(badChans == c, 1);
            reasons{idx} = [reasons{idx} ' + kurtosis z > 5'];
        end
    end
catch
end

[badChans, si] = sort(badChans);
reasons = reasons(si);
end

function interpChans = prompt_channel_interp(EEG, suggested)
fprintf('\n[INITREJ] Suggested bad channels: %s\n', vec2str(suggested));
if ~isempty(suggested)
    for k = 1:numel(suggested)
        ch = suggested(k);
        fprintf('   %d) %d (%s)\n', k, ch, safe_chan_label(EEG, ch));
    end
end
fprintf('\nType channel indices to INTERPOLATE (e.g., [1 2 17])\n');
fprintf('Default = none (Press Enter or type []).\n');
resp = input('Channels to interpolate: ', 's');
resp = strtrim(resp);

if isempty(resp) || strcmp(resp, '[]')
    interpChans = [];
    return;
end

interpChans = str2num(resp);
if ~isnumeric(interpChans)
    error('Manual interpolation must be numeric vector like [1 2 3] or [].');
end
interpChans = unique(interpChans(:))';
end

function make_initrej_plots(LOGS, subjid, EEG, metrics, badChans)
ensure_dir(LOGS);

save_hist(LOGS, subjid, metrics.chan_std, 'chan_std');
save_hist(LOGS, subjid, metrics.chan_rms, 'chan_rms');

save_bar(LOGS, subjid, metrics.chan_std, badChans, 'chan_std');
save_bar(LOGS, subjid, metrics.chan_rms, badChans, 'chan_rms');

if has_chanlocs(EEG)
    save_topo_metric(LOGS, subjid, EEG, metrics.chan_std, 'STD');
    save_topo_metric(LOGS, subjid, EEG, metrics.chan_rms, 'RMS');
end

save_channel_psd_overview(LOGS, subjid, EEG);

if ~isempty(badChans)
    save_channel_psd_badchans(LOGS, subjid, EEG, badChans);
end

% Label index map
fid = fopen(fullfile(LOGS, sprintf('sub-%03d_chalabels.txt', subjid)), 'w');
if fid > 0
    for k = 1:EEG.nbchan
        fprintf(fid, '%d\t%s\n', k, safe_chan_label(EEG, k));
    end
    fclose(fid);
end
end

function save_hist(LOGS, subjid, v, name)
h = figure('Visible', 'off');
histogram(v);
title(sprintf('sub-%03d %s histogram', subjid, name), 'Interpreter', 'none');
xlabel(name); ylabel('Count');
saveas(h, fullfile(LOGS, sprintf('sub-%03d_initrej_hist_%s.png', subjid, name)));
close(h);
end

function save_bar(LOGS, subjid, v, badChans, name)
h = figure('Visible', 'off');
bar(v);
title(sprintf('sub-%03d %s (suggested marked)', subjid, name), 'Interpreter', 'none');
xlabel('Channel index'); ylabel(name);
hold on;
for c = badChans(:)'
    xline(c, '--');
end
hold off;
saveas(h, fullfile(LOGS, sprintf('sub-%03d)initrej_bar_%s.png', subjid, name)));
close(h);
end

function save_topo_metric(LOGS, subjid, EEG, v, label)
h = figure('Visible', 'off');
topoplot(v, EEG.chanlocs, 'electrodes', 'on');
title(sprintf('sub-%03d %s topoplot', subjid, label), 'Interpreter', 'none');
saveas(h, fullfile(LOGS, sprintf('sub-%03d_initrej_topo_%s.png', subjid, lower(label))));
close(h);
end

function save_channel_psd_overview(LOGS, subjid, EEG)
fs = EEG.srate;
data = double(EEG.data);
nChan = size(data, 1);

win = round(fs * 2);
nover = round(win * 0.5);
nfft = max(2^nextpow2(win), win);

Pxx = zeros(nfft/2+1, nChan);
for ch = 1:nChan
    [pxx, f] = pwelch(data(ch, :), win, nover, nfft, fs);
    Pxx(:, ch) = 10*log10(pxx);
end

med = median(Pxx, 2);
q1 = prctile(Pxx, 25, 2);
q3 = prctile(Pxx, 75, 2);

h = figure('Visible', 'off');
plot(f, med); hold on;
plot(f, q1, ':'); plot(f, q3, ':');
xlim([0 80]);
xlabel('Hz'); ylabel('Power (dB)');
title(sprintf('sub-%03d channel PSD (median + IQR)', subjid), 'Interpreter', 'none');
legend({'Median', '25%', '75%'}, 'Location', 'northeast');
saveas(h, fullfile(LOGS, sprintf('sub-%03d_initrej_psd_overview.png', subjid)));
close(h);
end

function save_channel_psd_badchans(LOGS, subjid, EEG, badChans)
fs = EEG.srate;
data = double(EEG.data);

win = round(fs * 2);
nover = round(win * 0.5);
nfft = max(2^nextpow2(win), win);

h = figure('Visible', 'off'); hold on;
for ch = badChans(:)'
    [pxx, f] = pwelch(data(ch, :), win, nover, nfft, fs);
    plot(f, 10*log10(pxx));
end
xlim([0 80]);
xlabel('Hz'); ylabel('Power (dB)');
title(sprintf('sub-%03d PSD overlay (suggested channels %s)', subjid, vec2str(badChans)), 'Interpreter', 'none');
saveas(h, fullfile(LOGS, sprintf('sub-%03d_initrej_psd_badchans.png', subjid)));
close(h);
end

%% -------------------
% Channel PSD metrics
% --------------------
function chanPSD = compute_channel_psd_metrics(EEG)
fs = EEG.srate;
data = double(EEG.data);
nChan = size(data, 1);

win = round(fs * 2);
nover = round(win * 0.5);
nfft = max(2^nextpow2(win), win);

chanPSD = struct();
chanPSD.line_ratio  = zeros(nChan, 1);
chanPSD.hf_ratio    = zeros(nChan, 1);
chanPSD.drift_ratio = zeros(nChan, 1);
chanPSD.alpha_ratio = zeros(nChan, 1);

for ch = 1:nChan
    [pxx, f] = pwelch(data(ch, :), win, nover, nfft, fs);

    line = bp_psd(f, pxx, [59 61]);
    lo = bp_psd(f, pxx, [55 59]);
    hi = bp_psd(f, pxx, [61 65]);
    chanPSD.line_ratio(ch) = line / max(lo + hi, eps);

    hf = bp_psd(f, pxx, [20 40]);
    lf = bp_psd(f, pxx, [1 12]);
    chanPSD.hf_ratio(ch) = hf / max(lf, eps);

    drift = bp_psd(f, pxx, [1 2]);
    chanPSD.drift_ratio(ch) = drift / max(lf, eps);

    a = bp_psd(f, pxx, [8 12]);
    allp = bp_psd(f, pxx, [1 40]);
    chanPSD.alpha_ratio(ch) = a / max(allp, eps);
end
end

function p = bp_psd(f, pxx, band)
m = (f >= band(1) & f < band(2));
if ~any(m), p = 0; return; end
p = trapz(f(m), pxx(m));
end

function write_channel_psd_csv(LOGS, subjid, EEG, chanPSD)
ensure_dir(LOGS);
csvPath = fullfile(LOGS, sprintf('sub-%03d_chan_psd_metrics.csv', subjid));
fid = fopen(csvPath, 'w');
if fid < 0
    warning('Could not write channel PSD CSV: %s', csvPath);
    return;
end

fprintf(fid, 'chan_idx,label,line_ratio,hf_ratio,drift_ratio,alpha_ratio\n');
for ch = 1:EEG.nbchan
    fprintf(fid, '%d,%s,%.6f,%.6f,%.6f,%.6f\n', ...
        ch, safe_chan_label(EEG, ch), ...
        chanPSD.line_ratio(ch), chanPSD.hf_ratio(ch), chanPSD.drift_ratio(ch), chanPSD.alpha_ratio(ch));
end
fclose(fid);
end

function save_chan_psd_topos(LOGS, subjid, EEG, chanPSD)
if ~has_chanlocs(EEG); return; end
save_topo_metric(LOGS, subjid, EEG, chanPSD.line_ratio, 'LINE_RATIO');
save_topo_metric(LOGS, subjid, EEG, chanPSD.hf_ratio, 'HF_RATIO');
save_topo_metric(LOGS, subjid, EEG, chanPSD.drift_ratio, 'DRIFT_RATIO');
save_topo_metric(LOGS, subjid, EEG, chanPSD.alpha_ratio, 'ALPHA_RATIO');
end

%% ----------------------------------------------------
% ICA: training copy with optional bed segment removal
% -----------------------------------------------------
function [EEGtrain, segInfo] = make_ica_training_copy(EEG, cfg, logf)
EEGtrain = EEG;
segInfo = strict('rempoved', false, 'n_intervals', 0, 'pct_time', 0, 'intervals', []);

if ~isfield(cfg, 'preproc') || ~isfield(cfg.preproc, 'initrej') || ~isfield(cfg.preproc.initrej, 'badseg') || ~cfg.preproc.initrej.badseg.enabled
    logmsg(logf, '[ICA-TRAIN] badseg disabled; using full data.');
    return;
end

thr = cfg.preproc.initrej.badseg.threshold_uv;
logmsg(logf, '[ICA-TRAIN] Detecting bad segments (the = %.1f uV) for training copy.', thr);

x = double(EEG.data);
badSamp = any(abs(x) > thr, 1);
intervals = mask_to_intervals(badSamp);

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
EEGtrain = eeg_checkset(EEgtrain);

segInfo.removed = true;
logmsg(logf, '[ICA-TRAIN] Removed bad segments from training copy.');
end

function intervals = mask_to_intervals(mask)
masl - mask(:)';
if ~any(mask)
    intervals = [];
    return;
end
d = diff([0 mask 0]);
starts = find(d == 1);
ends   = find(d == -1) -1;
intervals = [starts(:) ends(:)];
end

function tf = prompt_yesno(prompt, defaultTF)
resp = input(prompt, 's');
resp = lower(strtrim(resp));
if isempty(resp)
    tf = defaultTF;
elseif any(strcmp(resp, {'y', 'yes'}))
    tf = true;
elseif any(strcmp(resp, {'n', 'no'}));
    tf = false;
else
    tf = defaultTF;
end
end

%% ----------------
% Event Validation
% -----------------
function validate_events_before_epoch(EEG, wanted, logf)
if isstring(wanted)
    wanted = cellstr(wanted);
end
if isempty(EEG.event)
    error('No EEG.event present; cannot epoch.');
end

types = {EEG.event.type};
typeStr = strings(1, numel(types));
for k = 1:numel(types)
    typeStr(k) = string(types{k});
end
u = unique(typeStr);

logmsg(logf, '[EVENTS] Unique event types: %s', numel(u), strjoin(u, ', '));
for k = 1:numel(u)
    logmsg(logf, '  [EVENTS] %s: %d', u(k), sum(typeStr == u(k)));
end

wantedStr = string(wanted);
hit = intersect(u, wantedStr);

if isempty(hit)
    error('None of requested event_types found: %s', strjoin(wantedStr, ', '));
else
    for k = 1:numel(hit)
        logmsg(logf, '[EVENTS] Will epoch on "%s" (n = %d)', hit9k0, sum9typeStr == hit(k)));
    end
end
end

%% -------------------------------
% ICLabel: suggest rejection list
% --------------------------------
function [suggestICs, reasons] = iclabel_suggest_reject(EEG, thr)
% Uses EEG.etc.ic_classification.ICLabel.classifications
% Order: [Brain Muscle Eye Heart LineNoise ChannelNoise Other]

suggestICs = [];
reasons = {};

if ~isfield(EEG, 'etc') || ~isfield(EEG.etc, 'ic_classification') || ~isfield(EEG.etc.ic_classification, 'ICLabel') || ~isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
    return;
end

C = EEG.etc.ic_classification.ICLabel.classifications;
if isempty(C), return; end

for ic = 1:size(C, 1)
    pBrain  = C(ic, 1);
    pMus    = C(ic, 2);
    pEye    = C(ic, 3);
    pHeart  = C(ic, 4);
    pLine   = C(ic, 5);
    pChan   = C(ic, 6);
    pOther  = C(ic, 7);

    hits = {};

    if isfield(thr, 'eye') && pEye >= thr.eye
        hits{end+1} = sprintf('eye = %.2f >= %.2f', pEye, thr.eye);
    end
    if isfield(thr, 'muscle') && pMus >= thr.muscle
        hits{end+1} = sprintf('muscle = %.2f >= %.2f', pMus, thr.muscle);
    end
    if isfield(thr, 'heart') && pHeart >= thr.heart
        hits{end+1} = sprintf('heart = %.2f >= %.2f', pHeart, thr.heart);
    end
    if isfield(thr, 'line_noise') && pLine >= thr.line_noise
        hits{end+1} = sprintf('line = %.2f >= %.2f', pLine, thr.line_noise);
    end
    if isfield(thr, 'channel_noise') && PChan >= thr.channel_noise
        hits{end+1} = sprintf('chanNoise = %.2f >= %.2f', pChan, thr.channel_noise);
    end

    if ~isempty(hits)
        suggestICs(end+1) = ic;
        reasons{end+1} = sprintf('IC %d: %s (brain = %.2f other = %.2f', ...
            ic, strjoin(hits, ', '), pBrain, pOther);
    end
end
end

%% -------------
% IC QC Packets
% --------------
function save_ic_qc_packets(QC, subjid, EEG, icList)
if isempty(icList), return; end
ensure_dir(QC);

outDir = fullfile(QC, sprintf('sub-%03d_icqc', subjid));
ensure_dir(outDir);

fs = EEG.srate;

W = EEG.icaweights * EEG.icasphere;
X = double(EEG.data(EEG.icachansind, :));
act = W * X;

C = [];
if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classifications') && isfield(EEG.etc.ic_classifications, 'ICLabel') && isfield(EEG.etc.ic_classifications.ICLabel, 'classifications')
    C = EEG.etc.ic_classification.ICLabel.classifications;
end

for ic = icList(:)'
    h = figure('Visible', 'off');

    subplot(2, 2, 1);
    if isfield(EEG, 'icawinv') && has_chanlocs(EEG)
        topoplot(EEG.icawinv(:, ic), EEG.chanlocs, 'electrodes', 'on');
        title(sprintf('IC %d scalp map', ic), 'Interpreter', 'none');
    end

    subplot(2, 2, 2);
    win = round(fs * 2);
    nover = round(win * 0.5);
    nfft = max(2^nextpow2(win), win);
    [pxx, f] = pwelch(act(ic, :), win, nover, nfft, fs);
    plot(f, 10*log10(pxx));
    xlim([0 80]);
    xlabel('Hz'); ylabel('Power (dB)');
    title('IC PSD', 'Interpreter', 'none');

    subplot(2, 2, [3 4]);
    nSamp = min(size(act, 2), round(fs * 10));
    plot((0:nSamp-1) / fs, act(ic, 1:nSamp));
    xlabel('Time (s)'); ylabel('a.u.');
    title('IC activation (first 10s)', 'Interpreter', 'none');

    if ~isempty(C) && size(C, 1) >= ic
        p = C(ic, :);
        ttl = sprintf('sub-%03d IC %d | B%.2f M%.2f E%.2f H%.2f L%.2f C%.2f O%.2f', ...
            subjid, ic, p(1), p(2), p(3), p(4), p(5), p(6), p(7));
    else
        ttl = sprintf('sub-%03d IC %d', subjid, ic);
    end
    sgtitle(ttl, 'Interpreter', 'none');

    saveas(h, fullfile(outDir, sprintf('sub-$03d_ic%03d_qc.png', subjid, ic)));
    close(h);
end
end

function removedICs = prompt_ic_reject(suggestICs)
fprintf('\n[ICREJ] ICLabel suggested ICs: %s\n', vec2str(suggestICs));
fprintf('Review QC figs in QC/sub-XXX_icqc/ before deciding.\n');
fprintf('Type IC indices to REMOVE (e.g., [1 3 7]).\n');
fprintf('Default = remove none (press Enter or type []).\n');

resp = input('ICs to remove: ', 's');
resp = strtrim(resp);

if isempty(resp) || strcmp(resp, '[]')
    removedICs = [];
    return;
end

removedICs = srt2num(resp);
if isempty(removedICs)
    removedICs = [];
else
    removedICs = unique(removedICs(:))';
end
end

%% --------------
% IC PSD metrics
% ---------------
function icMetrics = compute_ic_psd_metrics(EEG, icList)
if isempty(icList)
    icMetrics = strict('ic', {}, 'peak_hz', {}, 'bp', {}, 'hf_ratio', {}, 'line_ratio', {});
    return;
end

fs = EEG.srate;

W = EEG.icaweights * EEG.icasphere;
X = double(EEG.data(EEG.icachansind, :));
act = W * X;

win = round(fs * 2);
nover = round(win * 0.5);
nfft = max(2^nextpow2(win), win);

bands = struct( ...
    'delta', [1 4], ...
    'theta', [4 8], ...
    'alpha', [8 12], ...
    'beta', [13 30], ...
    'gamma', [30 45]);

icMetrics = repmat(struct('ic', NaN, 'peak_hz', NaN, 'bp', NaN, 'hf_ratio', NaN, 'line_ratio', NaN), numel(icList), 1);

for k = 1:numel(icList)
    ic = icList(k);

    [pxx, f] = pwelch(act(ic, :), win, nover, nfft, fs);

    bandMask = (f >= 0.5 & f <= 40);
    pBand = pxx(bandMask);
    fband = f(bandMask);
    [~, idx] = max(pband);
    peakHz = fband(idx);

    bp = struct();
    fn = fieldnames(bands);
    for j = 1:numel(fn)
        b = bands.(fn{j});
        bp.(fn{j}) = bandpower_from_psd(f, pxx, b);
    end

    hf = bandpower_from_psd(f, pxx, [20 40]);
    lf = bandpower_from_psd(f, pxx, [1 12]);
    hf_ratio = hf / max(lf, eps);

    line = bandpower_from_psd(f, pxx, [59 61]);
    lo = bandpower_from_psd(f, pxx, [55 59]);
    hi = bandpower_from_psd(f, pxx, [61 65]);
    line_ratio = line / max(lo + hi, eps);

    icMetrics(k).ic = ic;
    icMetrics(k).peak_jz = peakHz;
    icMetrics(k).bp = bp;
    icMetrics(k).hf_ratio = hf_ratio;
    icMetrics(k).line_ratio = line_ratio;
end
end

function p = bandpower_from_psd(f, pxx, band)
m = (f >= band(1) & f < band(2));
if ~any(m), p = 0; return; end
p = trapz(f(m), pxx(m));
end

function write_ic_metrics_csv(LOGS, subjid, icMetrics)
ensure_dir(LOGS);
csvPath = fullfile(LOGS, sprintf('sub-%03d_ic_psd_metrics.csv', subjid));
fid = fopen(csvPath, 'w');
if fid < 0
    warning('Could not write IC metrics CSV: %s', csvPath);
    return;
end

fprintf(fid, 'ic,peak_hz,delta,theta,alpha,beta,gamma,hf_ratio,line_ratio\n');
for k = 1:numel(icMetrics)
    bp = icMetrics(k).bp;
    fprintf(fid, '%d,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6f,%.6f\n', ...
        icMetrics(k).ic, icMetrics(k).peak_hz, ...
        bp.delta, bp.theta, bp.alpha, bp.beta, bp.gamma, ...
        icMetrics(k).hf_ratio, icMetrics(k).line_ratio);
end
fclose(fid);
end

function log_ic_metrics(logf, icMetrics)
if isempty(icMetrics), return; end
logmsg(logf, '[ICEMT] Per-IC PSD summary (suggested ICs):');
for k = 1:numel(icMetrics)
    bp = icMetrics(k).bp;
    logmsg(logf, '  IC %d | peak = %.2f Hz | d = %.2e t = %.2f a = %.2f b = %.2f g = %.2f | HF = %.2f | line = %.3f', ...
        icMetrics(k).ic, icMetrics(k).peak_hz, ...
        bp.delta, bp.theta, bp.alpha, bp.beta, bp.gamma, ...
        icMetrics(k).hf_ratio, icMetrics(k).line_ratio);
end
end