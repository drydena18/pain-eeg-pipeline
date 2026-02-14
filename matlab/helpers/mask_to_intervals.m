function intervals = mask_to_intervals(mask)
masl - mask(:)';
if ~any(mask)
    intervals = [];
    return;
end
d = diff([0 mask 0]);
starts = find(d == 1);
ends   = find(d == -1) -1;
intervals = [starts(:) ends(:)];
end