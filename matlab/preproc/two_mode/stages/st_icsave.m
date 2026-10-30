function [EEG, qc] = st_icsave(EEG, cfg, ~)
%ST_ICSAVE Mode 2b, first stage: save the post-rejection dataset to 06_ica
%   with preproc_core's tags (..._ica_iclabel.set), so the 07/08 files get
%   the same names downstream code already expects. The ICs were removed by
%   pp_mode2b from the saved decisions (EEG.etc.pp_ic_rejected).
rt  = cfg.rt;
ica = cfg.exp.preproc.ica;

tags = st_tags(EEG);
tags{end+1} = char(string(ica.tag));
if isfield(ica, 'iclabel') && isfield(ica.iclabel, 'enabled') && ica.iclabel.enabled && ...
        isfield(ica.iclabel, 'tag') && ~isempty(ica.iclabel.tag)
    tags{end+1} = char(string(ica.iclabel.tag));
end

rej = EEG.etc.pp_ic_rejected;
sug = EEG.etc.pp_ic_suggested;
if isempty(rej)
    logmsg(rt.logf, '[ICREJ] No ICs removed.');
else
    logmsg(rt.logf, '[ICREJ] Removed ICs: %s', vec2str(rej));
end
kept = setdiff(sug, rej);
if ~isempty(kept), logmsg(rt.logf, '[ICREJ] Suggested but NOT removed: %s', vec2str(kept)); end

EEG.etc.iclabel = struct();
EEG.etc.iclabel.suggestICs = sug;
EEG.etc.iclabel.removedICs = rej;

EEG = st_tags(EEG, tags);
save_stage(rt.ST.ICA, rt.P, rt.subjid, tags, EEG, rt.logf);
qc = struct('n_removed', numel(rej));
qc.saved = fullfile(rt.ST.ICA, char(rt.P.NAMING.fname(rt.subjid, tags, [])));
end