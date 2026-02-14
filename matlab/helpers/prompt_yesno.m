function tf = prompt_yesno(prompt, defaultTF)
resp = input(prompt, 's');
resp = lower(strtrim(resp));
if isempty(resp)
    tf = defaultTF;
elseif any(strcmp(resp, {'y', 'yes'}))
    tf = true;
elseif any(strcmp(resp, {'n', 'no'}));
    tf = false;
else
    tf = defaultTF;
end
end