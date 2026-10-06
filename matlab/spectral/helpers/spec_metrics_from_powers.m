function feat = spec_metrics_from_powers(ps, pf, pa, cog)
% SPEC_METRICS_FROM_POWERS  The 10 alpha feature families (incl. psi_cog)
% V 1.1.0
%
% V1.1.0: slow_alpha_frac removed; it equals (1 + sf_balance) / 2 exactly.
%
% Python twin: src_alpha_features._metrics_from_powers. NaN propagates.
% All inputs must have the same size ([nUnit x nTr], [1 x nTr], [nUnit x 1]).

e = 1e-12;
feat = struct();
feat.pow_slow_alpha  = ps;
feat.pow_fast_alpha  = pf;
feat.pow_alpha_total = pa;
feat.paf_cog_hz      = cog;
feat.sf_ratio        = ps ./ (pf + e);
feat.sf_logratio     = log(ps + e) - log(pf + e);
feat.sf_balance      = (ps - pf) ./ (ps + pf + e);
feat.rel_slow_alpha  = ps ./ (pa + e);
feat.rel_fast_alpha  = pf ./ (pa + e);
feat.psi_cog         = feat.sf_balance .* (cog - 10);
end
