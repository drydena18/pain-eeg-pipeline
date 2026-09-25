function labels = get_chan_labels(EEG, ctx)
% GET_CHAN_LABELS Return channel labels as a 1 x nbchan cellstr, or fail
% loudly.
%
%   labels = get_chan_labels(EEG);
%   labels = get_chan_labels(EEG, 'sub-003 stage09');

if nargin < 2 || isempty(ctx), ctx = ''; end
ctx = char(string(ctx));

if ~isfield(EEG, 'chanlocs')
    error('get_chan_labels:noField', '%s EEG has no chanlocs field.', ctx);
end

cl = EEG.chanlocs;

if ~isstruct(cl)
    error('get_chan_labels:notStruct', ...
        '%s EEG.chanlocs is class "%s" (size %s), not a struct array. nbchan = %d', ...
        ctx, class(cl), mat2str(size(cl)), EEG.nbchan);
end

if ~isfield(cl, 'labels')
    error('get_chan_labels:noLabels', ...
        '%s EEG.chanlocs has no "labels" field (fields: %s).', ...
        ctx, strjoin(fieldnames(cl)', ', '));
end

if numel(cl) ~= EEG.nbchan
    error('get_chan_labels:countMismatch', ...
        '%s numel(EEG.chanlocs) = %d but EEG.nbchan = %d.', ...
        ctx, numel(cl), EEG.nbchan);
end

labels = cellfun(@(x) char(string(x)), {cl.labels}, 'UniformOutput', false);

if any(cellfun(@isempty, labels))
    error('get_chan_labels:emptyLabel', ...
        '%s %d channel(s) have an empty label (indices: %s)', ...
        ctx, nnz(cellfun(@isempty, labels)), mat2str(find(cellfun(@isempty, labels))));
end

if numel(unique(upper(labels))) ~= numel(labels)
    error('get_chan_labels:duplicateLabel', ...
        '%s duplicate channel labels found.', ctx);
end
end