function log_ic_metrics(logf, icMetrics)
if isempty(icMetrics), return; end
logmsg(logf, '[ICEMT] Per-IC PSD summary (suggested ICs):');
for k = 1:numel(icMetrics)
    bp = icMetrics(k).bp;
    logmsg(logf, '  IC %d | peak = %.2f Hz | d = %.2e t = %.2f a = %.2f b = %.2f g = %.2f | HF = %.2f | line = %.3f', ...
        icMetrics(k).ic, icMetrics(k).peak_hz, ...
        bp.delta, bp.theta, bp.alpha, bp.beta, bp.gamma, ...
        icMetrics(k).hf_ratio, icMetrics(k).line_ratio);
end
end