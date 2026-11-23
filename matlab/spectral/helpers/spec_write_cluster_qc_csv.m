function spec_write_cluster_qc_csv(outPath, subjid, C, stimSide, minCh)
% SPEC_WRITE_CLUSTER_QC_CSV  Which electrodes made up each cluster, per subject
% V 1.0.0
%
% One row per cluster unit. Lets you check, across experiments with
% different caps, that every cluster was built from the same electrodes
% (and see which ones were missing where).
%
% Columns: subjid, cluster, kind, stim_side, n_requested, n_found,
%          min_channels, ok, requested, found, missing
%          (electrode lists are ';'-separated)

nC = numel(C);
T = table();
T.subjid       = repmat(int32(subjid), nC, 1);
T.cluster      = string({C.name})';
T.kind         = string({C.kind})';
T.stim_side    = repmat(string(stimSide), nC, 1);
T.n_requested  = arrayfun(@(c) numel(c.requested), C(:));
T.n_found      = arrayfun(@(c) numel(c.idx), C(:));
T.min_channels = repmat(double(minCh), nC, 1);
T.ok           = [C.ok]';
T.requested    = local_join({C.requested});
T.found        = local_join({C.found});
T.missing      = local_join({C.missing});
writetable(T, outPath);
end

function s = local_join(lists)
% ';'-joined electrode list per cluster, '' for an empty list
s = strings(numel(lists), 1);
for k = 1:numel(lists)
    v = lists{k};
    if ~isempty(v)
        s(k) = string(strjoin(reshape(cellstr(v), 1, []), ';'));
    end
end
end