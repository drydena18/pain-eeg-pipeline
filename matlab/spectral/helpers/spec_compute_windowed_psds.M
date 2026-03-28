function [f, prePxx, postPxx] = spec_compute_windowed_psds(EEG, psd, windows, logf)
% SPEC_COMPUTE_WINDOWED_PSDS  Pre- and post-stimulus Welch PSD for every trial
%
% Both windows share a common Welch window length and nfft so that the
% returned frequency vector f is identical for prePxx and postPxx.  This is
% required for the ERD computation in spec_compute_interaction_metrics.
%
% Inputs:
%   EEG     : EEGLAB struct (epoched; EEG.times in ms, t = 0 = stimulus onset)
%   psd     : psd config struct (window_sec, overlap_frac, nfft, fmin_hz, fmax_hz)
%   windows : struct with
%               pre_sec  [1x2]  e.g. [-1.0, -0.1] s relative to onset
%               post_sec [1x2]  e.g. [ 0.1,  0.8] s
%   logf    : (optional) MATLAB file handle for spec_logmsg
%
% Outputs:
%   f       : [1 x nFreq]
%   prePxx  : [nChan x nFreq x nTrials]  (NaN where a trial failed)
%   postPxx : [nChan x nFreq x nTrials]

if nargin < 4, logf = 1; end

% ---------------------------------------------------------------
% Build logical sample masks from EEG.times (ms -> s)
% ---------------------------------------------------------------
timesSec = double(EEG.times(:))' / 1000;

idxPre  = (timesSec >= windows.pre_sec(1))  & (timesSec <= windows.pre_sec(2));
idxPost = (timesSec >= windows.post_sec(1)) & (timesSec <= windows.post_sec(2));

nPre  = nnz(idxPre);
nPost = nnz(idxPost);

if nPre < 8
    error('spec_compute_windowed_psds:ShortPre', ...
        'Pre-stim window yields only %d samples (minimum 8). ' ...
        'Check cfg.spectral.windows.pre_sec and epoch length.', nPre);
end
if nPost < 8
    error('spec_compute_windowed_psds:ShortPost', ...
        'Post-stim window yields only %d samples (minimum 8). ' ...
        'Check cfg.spectral.windows.post_sec and epoch length.', nPost);
end

% ---------------------------------------------------------------
% Derive a COMMON Welch window / nfft so both windows share f
% ---------------------------------------------------------------
fs          = EEG.srate;
winS        = psd.window_sec;
if isempty(winS) || winS <= 0, winS = 2.0; end
win_desired = max(8, round(winS * fs));

% Must not exceed the shorter of the two windows
win_common = min(win_desired, min(nPre, nPost));
win_common = max(win_common, 8);

nfft_cfg = psd.nfft;
if isempty(nfft_cfg) || nfft_cfg <= 0
    nfft_common = max(2^nextpow2(win_common), win_common);
else
    nfft_common = max(nfft_cfg, win_common);
end

% Build a psd config override so spec_compute_psd_epoch_window uses the
% same win / nfft for both calls
psd_common              = psd;
psd_common.window_sec   = win_common / fs;  % ensures round(win_common/fs*fs)=win_common
psd_common.nfft         = nfft_common;

spec_logmsg(logf, '[WINDOWED_PSD] Common Welch: win=%d samp (%.3fs) nfft=%d | pre=%d post=%d samp', ...
    win_common, win_common/fs, nfft_common, nPre, nPost);

fmin = psd.fmin_hz;
fmax = psd.fmax_hz;

nChan = EEG.nbchan;
nTr   = EEG.trials;

% Pre-allocate after the first successful trial
f       = [];
prePxx  = [];
postPxx = [];

for t = 1:nTr
    try
        [fPre,  pPre ] = spec_compute_psd_epoch_window(EEG, psd_common, t, idxPre,  fmin, fmax);
        [~,     pPost] = spec_compute_psd_epoch_window(EEG, psd_common, t, idxPost, fmin, fmax);
    catch ME
        spec_logmsg(logf, '[WINDOWED_PSD][WARN] Trial %d failed: %s', t, ME.message);
        continue;
    end

    if isempty(f)
        f       = fPre;      % [1 x nFreq]
        nF      = numel(f);
        prePxx  = nan(nChan, nF, nTr);
        postPxx = nan(nChan, nF, nTr);
    end

    prePxx(:, :, t)  = pPre;
    postPxx(:, :, t) = pPost;
end

if isempty(f)
    error('spec_compute_windowed_psds:NoTrials', ...
        'No trials could be processed. Check epoch length vs window settings.');
end
end