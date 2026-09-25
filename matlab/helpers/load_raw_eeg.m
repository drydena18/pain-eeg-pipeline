function EEG = load_raw_eeg(rawPath, logf)
% LOAD_RAW_EEG Format-aware raw loader for the preprocessing pipeline.
% V 1.0.0
%
%   .set            -> pop_loadset
%   .vhdr           -> pop_loadbv
%   .eeg / .vmrk    -> redirected to sibling .vhdr
%   .bdf / .edf     -> pop_biosig

if nargin < 2, logf = 1; end
rawPath = char(string(rawPath));
[d, n, x] = fileparts(rawPath);

switch lower(x)
    case '.set'
        EEG = pop_loadset('filename', [n x], 'filepath', d);

    case {'.vhdr', '.eeg', '.vmrk'}
        hdr = [n '.vhdr'];
        if ~isfile(fullfile(d, hdr))
            error('load_raw_eg:NoVhdr', ...
                'BrainVision header not found next to %s (expected %s).', rawPath, hdr);
        end
        if ~exist('pop_loadbv', 'file')
            error('load_raw_eeg:NoBvaIo', ...
            ['pop_loadbv not on path. Install the "bva-io" plugin via ', ...
            'EEGLAB > File > Manage EEGLAB extensions.']);
        end
            logmsg(logf, '[LOAD] BrainVision via pop_loadbv: %s', fullfile(d, hdr));
            EEG = pop_loadbv([d filesep], hdr);

    case {'.bdf', '.edf'}
        logmsg(logf, '[LOAD] BIOSIG via pop_biosig: %s', rawPath);
        EEG = pop_biosig(rawPath);

    otherwise
       error('load_raw_eeg:UnknownExt', 'Unsupported raw extension "%s" (%s).', x,rawPath);
end

    EEG.etc.raw_file = rawPath;
end