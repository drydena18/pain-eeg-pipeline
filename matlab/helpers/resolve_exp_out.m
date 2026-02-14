function expOut = resolve_exp_out(P, cfg)
% Prefer registry out_dirname, else cfg, else exp_id
expOut = "";
if isfield(P, 'EXP') && isfield(P.EXP, 'out_dirname')
    expOut = string(P.EXP.out_dirname);
elseif isfield(cfg, 'exp') && isfield(cfg.exp, 'out_dirname')
    expOut = string(cfg.exp.out_dirname);
elseif isfield(cfg, 'exp') && isfield(cfg.exp, 'id')
    expOut = string(cfg.exp.id);
else
    expOut = "experiment";
end
end