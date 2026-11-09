function T = pp_diag_chanlocs(pp)
%PP_DIAG_CHANLOCS Where did EEG.chanlocs disappear? Read-only.
%   For every stage file in sub-XXX/00_concat ... 06_ica_ready, reads the
%   saved header directly (no eeg_checkset) and reports nbchan, the number
%   of channel locations, how many have coordinates, and when the file was
%   written - files older than this run were reused from an earlier
%   preproc_core run. Also prints the relevant lines from the Mode 1 log.
%
%   preproc_2mode('exp06', 'chanlocs', [55 59 66 73 95])

stages = {'00_concat', '01_filter', '02_notch', '03_resample', '04_reref', ...
          '05_initrej', '06_ica', '06_ica_ready'};
rows = {};
for i = 1:numel(pp.subjects)
    subj = pp.subjects{i};
    subDir = fullfile(pp.paths.out, subj);
    firstBad = '';
    for k = 1:numel(stages)
        d = dir(fullfile(subDir, stages{k}, '*.set'));
        for j = 1:numel(d)
            f = fullfile(d(j).folder, d(j).name);
            [nb, nl, nxyz, note] = local_header(f);
            if isempty(firstBad) && ~isnan(nb) && nl ~= nb
                firstBad = stages{k};
                note = strtrim([note ' <-- first mismatch']);
            end
            rows(end+1, :) = {subj, stages{k}, d(j).name, ...
                char(datetime(d(j).datenum, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm')), ...
                nb, nl, nxyz, note}; %#ok<AGROW>
        end
    end

    logf = fullfile(subDir, 'LOGS', sprintf('%s_pp_mode1.log', subj));
    if isfile(logf)
        L = splitlines(string(fileread(logf)));
        L = L(contains(L, ["[SKIP]", "[CHANTSV]", "Chanlocs", "chanlocs", "ELP", "[CHANLABELS]", ...
                           "[INITREJ] Interpolating", "[CONCAT] Merged", "Session"]));
        fprintf('\n--- %s: Mode 1 log (channel-related lines) ---\n', subj);
        fprintf('%s\n', L(max(1, end - 40):end));
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'subject', 'stage', 'file', 'written', 'nbchan', 'n_chanlocs', 'n_with_xyz', 'note'});
fprintf('\n');
disp(T);
end

% -------------------------------------------------------------------------
function [nb, nl, nxyz, note] = local_header(f)
nb = NaN; nl = NaN; nxyz = NaN; note = '';
try
    vars = who('-file', f);
    if ismember('EEG', vars)
        S = load(f, '-mat', 'EEG'); H = S.EEG;
    else
        H = load(f, '-mat', 'nbchan', 'chanlocs');
    end
    nb = double(H.nbchan);
    if isstruct(H.chanlocs)
        nl = numel(H.chanlocs);
        if isfield(H.chanlocs, 'X')
            nxyz = sum(arrayfun(@(c) ~isempty(c.X) && ~isnan(c.X), H.chanlocs));
        else
            nxyz = 0;
        end
    else
        nl = 0; nxyz = 0;
        note = sprintf('chanlocs is %s', class(H.chanlocs));
    end
catch ME
    note = ['unreadable: ' ME.message];
end
end
