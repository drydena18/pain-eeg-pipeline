function write_ic_metrics_csv(LOGS, subjid, icMetrics)
ensure_dir(LOGS);
csvPath = fullfile(LOGS, sprintf('sub-%03d_ic_psd_metrics.csv', subjid));
fid = fopen(csvPath, 'w');
if fid < 0
    warning('Could not write IC metrics CSV: %s', csvPath);
    return;
end

fprintf(fid, 'ic,peak_hz,delta,theta,alpha,beta,gamma,hf_ratio,line_ratio\n');
for k = 1:numel(icMetrics)
    bp = icMetrics(k).bp;
    fprintf(fid, '%d,%.4f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6f,%.6f\n', ...
        icMetrics(k).ic, icMetrics(k).peak_hz, ...
        bp.delta, bp.theta, bp.alpha, bp.beta, bp.gamma, ...
        icMetrics(k).hf_ratio, icMetrics(k).line_ratio);
end
fclose(fid);
end