function cfg = load_cfg(json_path)
% LOAD_CFG Load JSON configuration into a MATLAB struct

    if ~exist(json_path, 'file')
        error('Config file not found: %s.', json_path)
    end
    txt = fileread(json_path);
    cfg = jsondecode(txt);

end