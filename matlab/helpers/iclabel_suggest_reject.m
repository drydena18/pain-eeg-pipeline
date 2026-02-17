function [suggestICs, reasons] = iclabel_suggest_reject(EEG, thr)
% Uses EEG.etc.ic_classification.ICLabel.classifications
% Order: [Brain Muscle Eye Heart LineNoise ChannelNoise Other]

suggestICs = [];
reasons = {};

if ~isfield(EEG, 'etc') || ~isfield(EEG.etc, 'ic_classification') || ~isfield(EEG.etc.ic_classification, 'ICLabel') || ~isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
    return;
end

C = EEG.etc.ic_classification.ICLabel.classifications;
if isempty(C), return; end

for ic = 1:size(C, 1)
    pBrain  = C(ic, 1);
    pMus    = C(ic, 2);
    pEye    = C(ic, 3);
    pHeart  = C(ic, 4);
    pLine   = C(ic, 5);
    pChan   = C(ic, 6);
    pOther  = C(ic, 7);

    hits = {};

    if isfield(thr, 'eye') && pEye >= thr.eye
        hits{end+1} = sprintf('eye = %.2f >= %.2f', pEye, thr.eye);
    end
    if isfield(thr, 'muscle') && pMus >= thr.muscle
        hits{end+1} = sprintf('muscle = %.2f >= %.2f', pMus, thr.muscle);
    end
    if isfield(thr, 'heart') && pHeart >= thr.heart
        hits{end+1} = sprintf('heart = %.2f >= %.2f', pHeart, thr.heart);
    end
    if isfield(thr, 'line_noise') && pLine >= thr.line_noise
        hits{end+1} = sprintf('line = %.2f >= %.2f', pLine, thr.line_noise);
    end
    if isfield(thr, 'channel_noise') && pChan >= thr.channel_noise
        hits{end+1} = sprintf('chanNoise = %.2f >= %.2f', pChan, thr.channel_noise);
    end

    if ~isempty(hits)
        suggestICs(end+1) = ic;
        reasons{end+1} = sprintf('IC %d: %s (brain = %.2f other = %.2f', ...
            ic, strjoin(hits, ', '), pBrain, pOther);
    end
end
end