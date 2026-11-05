function pp = pp_hooks(pp)
%PP_HOOKS Wires the two-mode runner to the CNED preprocessing pipeline.
%   Each hook mirrors one block of preproc_core (V2.2.1) and calls the same
%   helpers, stage folders and file names (stages/st_*.m). Shared by every
%   experiment; experiment settings come from cfg.exp (the JSON).

pp.runtime_fn = @st_runtime;          % subject id, 00-09 folders, LOGS/QC, log file

% ---- Mode 1: stages 00-05, clean copy, rank-aware ICA -----------------------
pp.mode1.load_fn = @st_load;          % 00: concat / load, BIDS labels, montage, ELP, label guard
pp.mode1.stages = struct( ...
    'name', {'filter', 'notch', 'resample', 'reref', 'initrej'}, ...   % preproc_core order
    'fn',   {@st_filter, @st_notch, @st_resample, @st_reref, @st_initrej});
pp.mode1.clean_copy_fn = @st_clean_copy;

% Unattended replacements for preproc_core's keyboard prompts. These are
% methodological choices: they are in the Mode 1 config hash and are
% logged in EEG.etc.initrej / EEG.etc.pipeline and the ICA-ready marker.
%   interp_policy  'suggested' : interpolate every channel suggest_bad_channels flags
%                  'none'      : interpolate nothing (prompt_channel_interp's Enter default)
%   badseg_policy  'remove'    : drop >threshold segments from the ICA training copy
%                  'keep'      : train on all data (the prompt_yesno default)
pp.mode1.interp_policy = 'suggested';
pp.mode1.badseg_policy = 'keep';

% ---- Mode 2a: ICLabel + QC packets, IC PSD metrics, topo grid ---------------
pp.mode2a.prep_fn = @st_icqc;

% ---- Review: component scroll -> prompt_ic_reject -> confirm scroll ----------
pp.review.ui_fn = @st_review_ui;

% ---- Mode 2b: 06_ica -> 07_epoch -> 08_base -> 09_hilbert ----------------------
pp.mode2b.stages = struct( ...
    'name', {'icasave', 'epoch', 'baseline', 'hilbert'}, ...
    'fn',   {@st_icsave, @st_epoch, @st_baseline, @st_hilbert});
end