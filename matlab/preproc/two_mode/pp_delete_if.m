function pp_delete_if(file)
%PP_DELETE_IF Delete a file if it exists.
if isfile(file), delete(file); end
end