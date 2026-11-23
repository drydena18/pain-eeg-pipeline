function [S, labels] = spec_append_units(A, labelsA, B, labelsB)
% SPEC_APPEND_UNITS  Append extra units (e.g. clusters) below channel rows
% V 1.0.0
%
% For every field of A whose first dimension is numel(labelsA), the rows
% of B's same-named field are appended underneath. A field missing from B
% (or with the wrong number of rows or columns) gets NaN rows, so the
% CSV writers always see a rectangular [nA + nB x n] array. Fields only in
% B are dropped, which keeps the CSV columns identical to a channel-only
% run.
%
% Inputs:
%   A, labelsA : channel features ([nA x n] fields) and labels
%   B, labelsB : cluster features ([nB x n] fields) and names
%
% Outputs:
%   S      : struct, fields [nA + nB x n]
%   labels : column cellstr of nA + nB labels

la = local_cellstr(labelsA);
lb = local_cellstr(labelsB);
nA = numel(la);
nB = numel(lb);
labels = [la; lb];

S = struct();
fn = fieldnames(A);
for k = 1:numel(fn)
    a = A.(fn{k});
    if ~(isnumeric(a) || islogical(a)) || size(a, 1) ~= nA
        continue;
    end
    a = double(a);
    b = nan(nB, size(a, 2));
    if isstruct(B) && isfield(B, fn{k})
        v = double(B.(fn{k}));
        if isequal(size(v), [nB, size(a, 2)])
            b = v;
        end
    end
    S.(fn{k}) = [a; b];
end
end

function out = local_cellstr(x)
if ischar(x)
    out = {x};
elseif iscell(x)
    out = cellfun(@(v) char(v), x(:), 'UniformOutput', false);
else
    out = cellstr(x);
    out = out(:);
end
end
