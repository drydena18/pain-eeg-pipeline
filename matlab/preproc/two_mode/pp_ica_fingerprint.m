function fp = pp_ica_fingerprint(W, S)
%PP_ICA_FINGERPRINT Identifier for an ICA decomposition that survives
%   EEGLAB's own rescaling. eeg_checkset (called by pop_saveset/pop_loadset)
%   multiplies each row of icaweights by its component's RMS when
%   option_scaleicarms is on (the default), so the raw weights change on
%   save. The unmixing rows (W*S), each scaled to unit length, do not.
%
%   fp = pp_ica_fingerprint(EEG)   or   pp_ica_fingerprint(icaweights, icasphere)
if nargin == 1
    S = W.icasphere;
    W = W.icaweights;
end
U = double(W) * double(S);
U = U ./ sqrt(sum(U .^ 2, 2));
fp = pp_hash([size(U)'; round(U(:) * 1e6)]);
end
