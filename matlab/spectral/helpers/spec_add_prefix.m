function out = spec_add_prefix(s, prefix)
% SPEC_ADD_PREFIX Prefix every field name in a struct
% V 1.0.0
%
% Returns a new struct with each field of s renamed to [prefix fieldname].
% Used to tag whole-epoch / pre-stim / post-stim feature structs so they
% can be merged into a single struct without name collisions
%
% Inputs:
%   s       : struct
%   prefix  : char/string, e.g., "whole_", "pre_", "post_"
% 
% Output:
%   out     : struct with renamed fields, same values, same order as s

out = struct();
fn = fieldnames(s);
prefix = char(prefix);
for i = 1:numel(fn)
    out.([prefix fn{i}]) = s.(fn{i});
end
end