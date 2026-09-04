function prompt_scroll_eeg(EEG, block, logf, mode, titleStr)
% PROMPT_SCROLL_EEG Open an EEGLAB scroll browser for manual visual QC.
%   block       : the specific scroll config struct for this call site, e.g.,
%                 cfg.preproc.initrej.scroll, cfg.preproc.ica.train_scroll,
%                 cfg.preproc.ica.scroll, cfg.preproc.ica.confirm.scroll.
%                 Fields: enabled (required), winlength_sec, spacing_uv
%                 (spacing_uv only used in 'channels' mode).
%   mode        : 'channels' -> continuous channel-data scroll
%                 'components' -> IC activation scroll (required EEG.icaweights)
%   titleStr    : window title / log label for this call site, so multiple
%                 scroll points (which may share a mode) are distinguishable.
%
% No-op (returns immediately) if block.enabled is false/missing, so this
% is safe to call unconditionally from preproc_core.m

if nargin < 5 || isempty(titleStr)
    titleStr = sptrinf('%s scroll', mode);
end

if isempty(block) || ~isfield(block, 'enabled') || ~logical(block.enabled)
    return;
end

switch lower(mode)
    case 'channels'
        icacomp = 1;
    case 'components'
        icacomp = 0;
        if ~isfield(EEG, 'icaweights') || isempty(EEG.icaweights)
            logmsg(logf, '[WARN] Component scroll (%s) requested but EEG.icaweights is empty; skipping', titleStr);
            return;
        end
    otherwise
        error('prompt_scroll_eeg:BadMode', 'mode must be `channels` or `components`');
end

winlength = 10;
if isfield(block, 'winlength_sec') && ~isempty(block.winlength_sec)
    winlength = block.winlength_sec;
end

logmsg(logf, '[SCROLL] Opening %s (winlength = %gs). Close the window when done reviewing, then press Enter to continue.', titleStr, winlength);

try
    args = {'winlength', winlength, 'title', titleStr};
    if strcmpi(mode, 'channels') && isfield(block, 'spacing_uv') && ~isempty(block.spacing_uv)
        args = [args, {'spacing', block.spacing_uv}];
    end
    pop_eegplot(EEG, icacomp, 1, 0, 0, args{:});
    input('Press Enter once done reviewing the scroll to continue: ', 's');
catch ME
    logmsg(logf, '[WARN] Scroll (%s) failed: %s', titleStr, ME.message);
end
end