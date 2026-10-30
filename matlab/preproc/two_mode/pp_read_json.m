function s = pp_read_json(file)
%PP_READ_JSON
s = jsondecode(fileread(file));
end