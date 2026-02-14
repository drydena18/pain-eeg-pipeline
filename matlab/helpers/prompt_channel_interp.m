function interpChans = prompt_channel_interp(EEG, suggested)
fprintf('\n[INITREJ] Suggested bad channels: %s\n', vec2str(suggested));
if ~isempty(suggested)
    for k = 1:numel(suggested)
        ch = suggested(k);
        fprintf('   %d) %d (%s)\n', k, ch, safe_chan_label(EEG, ch));
    end
end
fprintf('\nType channel indices to INTERPOLATE (e.g., [1 2 17])\n');
fprintf('Default = none (Press Enter or type []).\n');
resp = input('Channels to interpolate: ', 's');
resp = strtrim(resp);

if isempty(resp) || strcmp(resp, '[]')
    interpChans = [];
    return;
end

interpChans = str2num(resp);
if ~isnumeric(interpChans)
    error('Manual interpolation must be numeric vector like [1 2 3] or [].');
end
interpChans = unique(interpChans(:))';
end