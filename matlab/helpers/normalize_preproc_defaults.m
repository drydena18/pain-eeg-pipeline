function Pp = normalize_preproc_defaults(Pp)

% ---- Stage 00: multi-session concatenation ----
Pp = ensureBlock(Pp, 'concat', false, 'concat');
Pp.concat = defaultField(Pp.concat, 'n_sessions',      1);
Pp.concat = defaultField(Pp.concat, 'session_pattern', '');

% ---- Stage 01: highpass/lowpass filter ----
Pp = ensureBlock(Pp, 'filter', true, 'fir');
Pp.filter = defaultField(Pp.filter, 'type',         'fir');
Pp.filter = defaultField(Pp.filter, 'highpass_hz',  0.5);
Pp.filter = defaultField(Pp.filter, 'lowpass_hz',   40);

% ---- Stage 02: notch (bandstop) ----
Pp = ensureBlock(Pp, 'notch', true, 'notch60');
Pp.notch = defaultField(Pp.notch, 'freq_hz', 60);
Pp.notch = defaultField(Pp.notch, 'bw_hz',   2);

% ---- Stage 03: resample ----
Pp = ensureBlock(Pp, 'resample', false, 'rs500');
Pp.resample = defaultField(Pp.resample, 'target_hz', []);

% ---- Stage 04: re-reference ----
Pp = ensureBlock(Pp, 'reref', true, 'reref');
Pp.reref = defaultField(Pp.reref, 'mode',     'average');
Pp.reref = defaultField(Pp.reref, 'channels', []);

% ---- Stage 05: initial rejection (bad channels + bad segments) ----
Pp = ensureBlock(Pp, 'initrej', true, 'initrej');
Pp.initrej = defaultStruct(Pp.initrej, 'badchan');
Pp.initrej = defaultStruct(Pp.initrej, 'badseg');
Pp.initrej.badchan = defaultField(Pp.initrej.badchan, 'enabled', false);
Pp.initrej.badseg  = defaultField(Pp.initrej.badseg,  'enabled', false);

% ---- Stage 06: ICA ----
Pp = ensureBlock(Pp, 'ica', true, 'ica');
Pp.ica = defaultField(Pp.ica, 'method', 'runica');
Pp.ica = defaultStruct(Pp.ica, 'iclabel');
Pp.ica.iclabel = defaultField(Pp.ica.iclabel, 'enabled',     false);
Pp.ica.iclabel = defaultField(Pp.ica.iclabel, 'auto_reject', false);
Pp.ica.iclabel = defaultStruct(Pp.ica.iclabel, 'thresholds');

% ---- Stage 07: epoch ----
Pp = ensureBlock(Pp, 'epoch', true, 'epoch');
Pp.epoch = defaultField(Pp.epoch, 'events', {});
Pp.epoch = defaultField(Pp.epoch, 'tmin_sec',    -1.0);
Pp.epoch = defaultField(Pp.epoch, 'tmax_sec',    2.0);

% ---- Stage 08: baseline correction ----
Pp = ensureBlock(Pp, 'baseline', true, 'base');
Pp.baseline = defaultField(Pp.baseline, 'window_sec', [-0.5 0]);

% Stage 09: Hilbert slow-alpha instantaneous phase
% Disabled by default. Enable in JSON: "hilbert": {"enabled": true}
% slow_hz should match cfg.spectral.alpha.slow_hz (both default to [8 10]).
% preproc_core reads cfg.spectral.alpha.slow_hz first if available, so this
% value acts only as a fallback when the spectral block is absent from the JSON.
Pp = ensureBlock(Pp, "hilbert", false, "hilbert");
Pp.hilbert = defaultField(Pp.hilbert, "slow_hz", [8 10]);

if Pp.epoch.enabled && isempty(Pp.epoch.events)
    error('preproc_default:EpochMissingEvents', ...
        'cfg.preproc.epoch.enabled is true, but epoch.events is empty.');
end
end