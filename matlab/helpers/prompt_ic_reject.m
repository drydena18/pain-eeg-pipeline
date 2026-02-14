function removedICs = prompt_ic_reject(suggestICs)
fprintf('\n[ICREJ] ICLabel suggested ICs: %s\n', vec2str(suggestICs));
fprintf('Review QC figs in QC/sub-XXX_icqc/ before deciding.\n');
fprintf('Type IC indices to REMOVE (e.g., [1 3 7]).\n');
fprintf('Default = remove none (press Enter or type []).\n');

resp = input('ICs to remove: ', 's');
resp = strtrim(resp);

if isempty(resp) || strcmp(resp, '[]')
    removedICs = [];
    return;
end

removedICs = srt2num(resp);
if isempty(removedICs)
    removedICs = [];
else
    removedICs = unique(removedICs(:))';
end
end