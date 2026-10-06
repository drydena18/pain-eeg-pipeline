function check_standard_chan_labels(EEG, logf, contextStr)
% CHECK_STANDARD_CHAN_LABELS  Stop if raw BioSemi / non-scalp labels survive import
% V 1.0.0
%
% Called by preproc_core after montage + channel-location lookup. Errors when
% any label still looks like a native BioSemi electrode (A1-A32, B1-B32) or an
% external/auxiliary channel (EXG*, Status, Erg*, GSR*, Resp*, Plet*, Temp*),
% because then:
%   - the .elp lookup cannot assign coordinates (interpolation, topoplots
%     and source localisation all need them), and
%   - non-scalp channels would enter the average reference.
% Fix: enable cfg.exp.montage (as in exp01) or provide a BIDS channels.tsv.
%
% Also logs (does not stop on) channels that have no 3-D coordinates, e.g.
% EOG channels in BrainVision recordings.

if nargin < 2, logf = 1; end
if nargin < 3, contextStr = ''; end

labs = {EEG.chanlocs.labels};
U = upper(strtrim(labs));

isAB  = ~cellfun(@isempty, regexp(U, '^(A|B)0*\d{1,2}$', 'once'));
isAux = ~cellfun(@isempty, regexp(U, '^(EXG\d*|STATUS|ERG\d*|GSR\d*|RESP|PLET|TEMP)$', 'once'));

bad = labs(isAB | isAux);
if ~isempty(bad)
    msg = sprintf(['%s %d channel(s) still carry raw BioSemi / auxiliary labels (%s). ', ...
        'Enable cfg.exp.montage (AB-select + CSV relabel, as in exp01) or supply a BIDS ', ...
        'channels.tsv before channel lookup.'], contextStr, numel(bad), strjoin(bad(1:min(10, end)), ', '));
    logmsg(logf, '[CHANLABELS][ERROR] %s', msg);
    error('check_standard_chan_labels:NonStandardLabels', '%s', msg);
end

noXYZ = false(1, numel(labs));
if isfield(EEG.chanlocs, 'X')
    noXYZ = arrayfun(@(c) isempty(c.X) || isnan(c.X), EEG.chanlocs);
end
if any(noXYZ)
    logmsg(logf, '[CHANLABELS][WARN] %s %d channel(s) have no coordinates: %s', ...
        contextStr, nnz(noXYZ), strjoin(labs(noXYZ), ', '));
end
end
