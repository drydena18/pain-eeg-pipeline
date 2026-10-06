function L = spec_lep_features(X, timesSec, fs, n2Win, p2Win)
% SPEC_LEP_FEATURES  N2 / P2 peak features per unit x trial
% V 1.0.0
%
% Python twin: src_lep._lep_arrays. Same windows, same rule, same names:
%   n2_amp / n2_lat_ms : most negative value in n2Win (first occurrence)
%   p2_amp / p2_lat_ms : most positive value in p2Win (first occurrence)
%   n2p2_amp           : p2_amp - n2_amp
%   n2_mean / p2_mean  : mean amplitude in each window
%
% Inputs:
%   X : [nUnit x nTime x nTr] (pass mean(X, 3) for trial-mean waveform)
%
% Output:
%   L : struct of [nUnit x nTr] arrays

timesSec = double(timesSec(:))';
nU = size(X, 1);
nTr = size(X, 3);
L = struct();
specs = {'n2', n2Win, 'min'; 'p2', p2Win, 'max'};
for r = 1:2
    nm = specs{r, 1};
    win = specs{r, 2};
    m = spec_window_mask(timesSec, win(1), win(2), fs);
    if ~any(m)
        L.([nm '_amp'])    = nan(nU, nTr);
        L.([nm '_lat_ms']) = nan(nU, nTr);
        L.([nm '_mean'])   = nan(nU, nTr);
        continue;
    end
    Xw = double(X(:, m, :));
    tw = timesSec(m);
    if strcmp(specs{r, 3}, 'min')
        [a, ix] = min(Xw, [], 2);
    else
        [a, ix] = max(Xw, [], 2);
    end
    L.([nm '_amp'])    = reshape(a, nU, nTr);
    L.([nm '_lat_ms']) = reshape(tw(ix) * 1000, nU, nTr);
    L.([nm '_mean'])   = reshape(mean(Xw, 2), nU, nTr);
end
L.n2p2_amp = L.p2_amp - L.n2_amp;
L = orderfields(L, {'n2_amp', 'n2_lat_ms', 'p2_amp', 'p2_lat_ms', 'n2p2_amp', 'n2_mean', 'p2_mean'});
end
