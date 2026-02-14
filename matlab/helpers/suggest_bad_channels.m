function [badChans, reasons, metrics] = suggest_bad_channels(EEG)
% Conservative automated suggestions using EEGLAB pop_rejchan prob + kurt
badChans = [];
reasons = {};

metrics = struct();
metrics.chan_rms = sqrt(mean(double(EEG.data).^2, 2));
metrics.chan_std = std(double(EEG.data), 0, 2);

try
    [~, badP] = pop_rejchan(EEG, 'threshold', 5, 'norm', 'on', 'measure', 'prob');
    for c = badP(:)'
        badChans(end+1) = c;
        reasons{end+1} = 'probability z > 5 (pop_rejchan)';
    end
catch
end

try
    [~, badK] = pop_rejchan(EEG, 'threshold', 5, 'norm', 'on', 'measure', 'kurt');
    for c = badK(:)'
        if ~ismember(c, badChans)
            badChans(end+1) = c;
            reasons{end+1} = 'kurtosis z > 5 (pop_rejchan)';
        else
            idx = find(badChans == c, 1);
            reasons{idx} = [reasons{idx} ' + kurtosis z > 5'];
        end
    end
catch
end

[badChans, si] = sort(badChans);
reasons = reasons(si);
end