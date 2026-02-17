function intervals = mask_to_intervals(mask)
mask = logical(mask);
mask = mask(:)';
d = diff([false mask false]);

starts = find(d == 1);
ends = find(d == -1) -1;

intervals = [starts(:) ends(:)];
end