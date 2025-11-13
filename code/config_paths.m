% Centralized paths for 26ByBiosemi EEG Pipeline
function P = config_paths()

% ---- Root ----
P.PROJ_ROOT = fullfile('/cifs','seminowicz','eegPainDatasets','CNED');

% ---- BIDS dataset root ----
P.BIDS = fullfile(P.PROJ_ROOT, '26ByBiosemi');

% ---- Core Inputs ----
P.PARTICIPANTS_TSV = fullfile(P.BIDS, 'participants.tsv');
P.DATASET_EVENTS_JSON = fullfile(P.BIDS, 'task-26ByBiosemi_events.json');
P.CSV_SINGLETRIAL = fullfile(P.PROJ_ROOT, 'data', 'CSV', 'participants_singletrial_26ByBiosemi.csv'); %% Need to downlaod

% ---- Derivatives ----
P.DERIV = fullfile(P.PROJ_ROOT,'da-analysis','26ByBiosemi','derivatives');
P.SESSION_MERGED = fullfile(P.DERIV, 'session_merged_data');
P.FILTER = fullfile(P.DERIV, 'filter');
P.NO_ICA = fullfile(P.DERIV, 'no_ica');
P.ICA = fullfile(P.DERIV, 'ica');
P.AFTER_ICA = fullfile(P.DERIV, 'after_ica');
P.REFILTER_30 = fullfile(P.DERIV, 'refilter_30');
P.REREFER = fullfile(P.DERIV, 'rerefer');

% ---- Resources ----
P.ELP_FILE = fullfile(P.PROJ_ROOT, 'resources', 'standard-10-5-cap385.elp'); %% Need to download

% Ensure derivative folders exist ----
outDirs = {P.SESSION_MERGED, P.FILTER, P.NO_ICA, P.ICA, P.AFTER_ICA, P.REFILTER_30, P.REREFER};
cellfun(@(d) ~exist(d, 'dir') && mkdir(d), outDirs);
end