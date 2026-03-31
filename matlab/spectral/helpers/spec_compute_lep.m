function lep = spec_compute_lep(EEG, lepCfg, outDir, subjid, logf)
% SPEC_COMPUTE_LEP  Extract Laser Evoked Potential (LEP) per trial and GA
% V 1.0.0
%
% Because stage 08_base is already baseline-corrected, the trial-level
% LEP is simply EEG.data windowed to the analysis epoch.  No additional
% baseline subtraction is applied here (it would double-correct).
%
% Saves three files under outDir:
%   sub-XXX_lep_trials.mat   : .data [nChan x nT x nTrials], .t_ms, .chan_labels, .fs
%   sub-XXX_lep_ga.mat       : .data [nChan x nT],            .t_ms, .chan_labels, .fs
%   sub-XXX_lep_peaks.csv    : per trial x channel N2 / P2 peak metrics
%
% Returns:
%   lep.ga_data      [nChan x nT]        grand-average waveform
%   lep.t_ms         [1 x nT]            time axis (ms, stimulus onset at 0)
%   lep.chan_labels   cell / string array
%   lep.peak_table   MATLAB table        (also written to CSV)
%
% Inputs:
%   EEG     : EEGLAB epoched struct (stage 08_base; EEG.times in ms)
%   lepCfg  : cfg.spectral.lep struct:
%               window_sec    [1x2]  analysis epoch, e.g. [-0.1, 1.0]
%               n2_window_sec [1x2]  N2 search,      e.g. [ 0.10, 0.30]
%               p2_window_sec [1x2]  P2 search,      e.g. [ 0.25, 0.45]
%   outDir  : folder to save .mat and .csv files
%   subjid  : integer subject ID
%   logf    : (optional) MATLAB file handle for spec_logmsg

if nargin < 5, logf = 1; end

winSec   = lepCfg.window_sec;     % e.g. [-0.1, 1.0] s
n2WinSec = lepCfg.n2_window_sec;  % e.g. [0.10, 0.30] s
p2WinSec = lepCfg.p2_window_sec;  % e.g. [0.25, 0.45] s

% ---------------------------------------------------------------
% Build sample index for the analysis window
% ---------------------------------------------------------------
timesSec = double(EEG.times(:))' / 1000;   % ms -> s
idxWin   = (timesSec >= winSec(1)) & (timesSec <= winSec(2));

if nnz(idxWin) < 2
    error('spec_compute_lep:EmptyWindow', ...
        'LEP window [%.3f, %.3f] s produces no samples in this epoch.', ...
        winSec(1), winSec(2));
end

t_ms       = EEG.times(idxWin);    % [1 x nT]
chanLabels = spec_get_chanlabels(EEG);
nChan      = EEG.nbchan;
nTr        = EEG.trials;

spec_logmsg(logf, '[LEP] window [%.3f %.3f] s -> %d samples, %d trials, %d chans', ...
    winSec(1), winSec(2), numel(t_ms), nTr, nChan);

% ---------------------------------------------------------------
% Windowed trial data  [nChan x nT x nTrials]  (double)
% ---------------------------------------------------------------
trialData = double(EEG.data(:, idxWin, :));

% Grand average across trials
gaData = mean(trialData, 3, 'omitnan');    % [nChan x nT]

% ---------------------------------------------------------------
% Save waveform .mat files
% ---------------------------------------------------------------
lepTrialStruct              = struct();
lepTrialStruct.data         = trialData;
lepTrialStruct.t_ms         = t_ms;
lepTrialStruct.chan_labels  = chanLabels;
lepTrialStruct.fs           = EEG.srate;
lepTrialStruct.subjid       = subjid;

lepGAStruct                 = struct();
lepGAStruct.data            = gaData;
lepGAStruct.t_ms            = t_ms;
lepGAStruct.chan_labels     = chanLabels;
lepGAStruct.fs              = EEG.srate;
lepGAStruct.subjid          = subjid;

trialMatPath = fullfile(outDir, sprintf('sub-%03d_lep_trials.mat', subjid));
gaMatPath    = fullfile(outDir, sprintf('sub-%03d_lep_ga.mat',     subjid));

save(trialMatPath, '-struct', 'lepTrialStruct');
save(gaMatPath,    '-struct', 'lepGAStruct');

spec_logmsg(logf, '[LEP] Saved trial .mat -> %s', trialMatPath);
spec_logmsg(logf, '[LEP] Saved GA    .mat -> %s', gaMatPath);

% ---------------------------------------------------------------
% Peak detection per trial x channel
%   N2 : most negative deflection inside n2_window_sec
%   P2 : most positive deflection inside p2_window_sec
% ---------------------------------------------------------------
idxN2 = (t_ms >= n2WinSec(1)*1000) & (t_ms <= n2WinSec(2)*1000);
idxP2 = (t_ms >= p2WinSec(1)*1000) & (t_ms <= p2WinSec(2)*1000);
t_n2  = t_ms(idxN2);
t_p2  = t_ms(idxP2);

if ~any(idxN2)
    spec_logmsg(logf, '[LEP][WARN] N2 search window [%.0f %.0f] ms empty; N2 peaks will be NaN.', ...
        n2WinSec(1)*1000, n2WinSec(2)*1000);
end
if ~any(idxP2)
    spec_logmsg(logf, '[LEP][WARN] P2 search window [%.0f %.0f] ms empty; P2 peaks will be NaN.', ...
        p2WinSec(1)*1000, p2WinSec(2)*1000);
end

nRows      = nChan * nTr;
subjCol    = repmat(subjid,             nRows, 1);
trialCol   = repelem((1:nTr)',          nChan);
chanIdxCol = repmat((1:nChan)',         nTr,   1);
chanLabCol = repmat(string(chanLabels(:)), nTr, 1);
n2Amp      = nan(nRows, 1);
n2Lat      = nan(nRows, 1);
p2Amp      = nan(nRows, 1);
p2Lat      = nan(nRows, 1);
n2p2Amp    = nan(nRows, 1);

for t = 1:nTr
    for ch = 1:nChan
        row = (t - 1)*nChan + ch;
        sig = squeeze(trialData(ch, :, t));   % [1 x nT]

        if any(idxN2)
            segN2       = sig(idxN2);
            [na, ni]    = min(segN2);
            n2Amp(row)  = na;
            n2Lat(row)  = t_n2(ni);
        end

        if any(idxP2)
            segP2       = sig(idxP2);
            [pa, pi_]   = max(segP2);
            p2Amp(row)  = pa;
            p2Lat(row)  = t_p2(pi_);
        end

        if any(idxN2) && any(idxP2)
            n2p2Amp(row) = p2Amp(row) - n2Amp(row);
        end
    end
end

peakTable = table( ...
    subjCol, trialCol, chanIdxCol, chanLabCol, ...
    n2Amp, n2Lat, p2Amp, p2Lat, n2p2Amp, ...
    'VariableNames', { ...
    'subjid', 'trial', 'chan_idx', 'chan_label', ...
    'n2_amp_uv', 'n2_lat_ms', 'p2_amp_uv', 'p2_lat_ms', 'n2p2_amp_uv'} ...
);

peakCSVPath = fullfile(outDir, sprintf('sub-%03d_lep_peaks.csv', subjid));
writetable(peakTable, peakCSVPath);
spec_logmsg(logf, '[LEP] Saved peak CSV -> %s', peakCSVPath);

% ---------------------------------------------------------------
% Pack return struct
% ---------------------------------------------------------------
lep.ga_data    = gaData;
lep.t_ms       = t_ms;
lep.chan_labels = chanLabels;
lep.peak_table  = peakTable;
end