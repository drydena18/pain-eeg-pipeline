function validate_events_before_epoch(EEG, wanted, logf)
if isstring(wanted)
    wanted = cellstr(wanted);
end
if isempty(EEG.event)
    error('No EEG.event present; cannot epoch.');
end

types = {EEG.event.type};
typeStr = strings(1, numel(types));
for k = 1:numel(types)
    typeStr(k) = string(types{k});
end
u = unique(typeStr);

logmsg(logf, '[EVENTS] Unique event types: %s', numel(u), strjoin(u, ', '));
for k = 1:numel(u)
    logmsg(logf, '  [EVENTS] %s: %d', u(k), sum(typeStr == u(k)));
end

wantedStr = string(wanted);
hit = intersect(u, wantedStr);

if isempty(hit)
    error('None of requested event_types found: %s', strjoin(wantedStr, ', '));
else
    for k = 1:numel(hit)
        logmsg(logf, '[EVENTS] Will epoch on "%s" (n = %d)', hit(k), sum(typeStr == hit(k)));
    end
end
end